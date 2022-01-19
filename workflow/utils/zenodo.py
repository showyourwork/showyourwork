"""
Main Zenodo interface.

"""
from . import exceptions, paths
from .logging import get_logger
from pathlib import Path
import requests
import os
import json
import subprocess


# Default evironment variable names
default_token_name = {
    "zenodo": "ZENODO_TOKEN",
    "zenodo_sandbox": "ZENODO_SANDBOX_TOKEN",
}

# Zenodo base URLs
zenodo_url = {"zenodo": "zenodo.org", "zenodo_sandbox": "sandbox.zenodo.org"}


def check_status(r):
    """
    Parse a requests return object and raise a custom exception
    for a >200-level status code.

    """
    if r.status_code > 204:
        data = r.json()
        for error in data.get("errors", []):
            data["message"] += " " + error["message"]
        raise exceptions.ZenodoError(data)
    return r


def get_access_token(token_name, error_if_missing=False):
    """
    Return the access token stored in the environment variable `token_name`.

    """
    access_token = os.getenv(token_name, None)
    if error_if_missing and access_token is None:
        raise exceptions.MissingZenodoAccessToken(token_name)
    return access_token


def _get_id_type(
    deposit_id,
    zenodo_url="zenodo.org",
    token_name="ZENODO_TOKEN",
    max_pages=100,
    results_per_page=10,
):
    """
    Determines whether a given Zenodo `id` corresponds to
    a concept id or a version id (and raises an error otherwise).

    """
    # Try to find a published record (no authentication needed)
    r = requests.get(f"https://{zenodo_url}/api/records/{deposit_id}")
    data = r.json()

    if r.status_code > 204:

        if "PID is not registered" in data.get("message", ""):

            # Get the access token
            access_token = get_access_token(token_name, error_if_missing=True)

            # There's no public record; let's search for a draft
            # deposit with this concept or version ID (authentication
            # needed)
            for page in range(1, max_pages + 1):
                r = check_status(
                    requests.get(
                        f"https://{zenodo_url}/api/deposit/depositions",
                        params={
                            "q": "",
                            "access_token": access_token,
                            "status": "draft",
                            "size": results_per_page,
                            "page": page,
                            "sort": "mostrecent",
                            "all_versions": 1,
                        },
                    )
                )
                data = r.json()
                for deposit in data:
                    if int(deposit["conceptrecid"]) == deposit_id:

                        # We found a deposit draft with this concept id
                        return "concept"

                    elif int(deposit["id"]) == deposit_id:

                        # We found a deposit draft with this version id
                        return "version"

                if len(data) == 0:

                    # We reached the end and found nothing
                    raise exceptions.ZenodoRecordNotFound(deposit_id)

        else:

            # Something unexpected happened...
            raise exceptions.ZenodoError(data)

    else:

        # This is a public record
        if int(deposit_id) == int(data["conceptrecid"]):
            return "concept"
        elif int(deposit_id) == int(data["id"]):
            return "version"
        else:
            raise exceptions.ZenodoRecordNotFound(deposit_id)


def get_id_type(deposit_id, zenodo_url="zenodo.org", **kwargs):
    """
    Determines whether a given Zenodo `id` corresponds to
    a concept id or a version id (and raises an error otherwise).

    Caches the result locally.

    """
    cache_file = paths.zenodo / f"{deposit_id}.{zenodo_url}"
    if cache_file.exists():

        with open(cache_file, "r") as f:
            id_type = f.readline().replace("\n", "")

    else:

        id_type = _get_id_type(deposit_id, zenodo_url=zenodo_url, **kwargs)
        with open(cache_file, "w") as f:
            print(id_type, file=f)

    return id_type


def upload_simulation(
    file_name,
    concept_id,
    deposit_title,
    deposit_description,
    deposit_creators,
    zenodo_url="zenodo.org",
    token_name="ZENODO_TOKEN",
    file_path=".",
    shell_cmd="",
    repo_url="",
    max_pages=100,
    results_per_page=10,
):
    # Set up the logger
    logger = get_logger()

    # Retrieve the access token
    access_token = get_access_token(token_name, error_if_missing=True)

    # Search for a public record with the correct concept id
    r = requests.get(f"https://{zenodo_url}/api/records/{concept_id}")

    if r.status_code > 204 and "PID is not registered" in r.json().get("message", ""):

        # There's no public record; let's search for a draft
        # deposit with this concept id (authentication needed)
        try:
            for page in range(1, max_pages + 1):
                r = check_status(
                    requests.get(
                        f"https://{zenodo_url}/api/deposit/depositions",
                        params={
                            "q": "",
                            "access_token": access_token,
                            "status": "draft",
                            "size": results_per_page,
                            "page": page,
                            "sort": "mostrecent",
                            "all_versions": 1,
                        },
                    )
                )
                data = r.json()
                for deposit in data:
                    if int(deposit["conceptrecid"]) == concept_id:
                        draft_id = int(deposit["id"])
                        raise StopIteration

                if len(data) == 0:
                    # We reached the end and found nothing
                    raise exceptions.ZenodoRecordNotFound(concept_id, "concept")

        except StopIteration:

            # We're good!
            pass

        else:

            # We reached the end and found nothing
            raise exceptions.ZenodoRecordNotFound(concept_id, "concept")

    elif r.status_code <= 204:

        # There's a public record; let's search for an existing draft
        data = r.json()
        latest_draft = data["links"].get("latest_draft", None)
        if latest_draft:

            # A draft exists; get its id
            draft_id = int(latest_draft.split("/")[-1])

        else:

            # There's no draft, so we need to create one
            latest_id = data["id"]
            r = check_status(
                requests.post(
                    f"https://{zenodo_url}/api/deposit/depositions/{latest_id}/actions/newversion",
                    params={"access_token": access_token},
                )
            )
            draft_id = int(r.json()["links"]["latest_draft"].split("/")[-1])

    else:

        raise exceptions.ZenodoRecordNotFound(concept_id, "concept")

    # Get the bucket url for this draft so we can upload files
    r = check_status(
        requests.get(
            f"https://{zenodo_url}/api/deposit/depositions/{draft_id}",
            params={"access_token": access_token},
        )
    )
    bucket_url = r.json()["links"]["bucket"]

    # Grab the latest version doi if we need to fall back to it below
    try:
        latest_id = r.json()["links"]["latest_html"].split("/")[-1]
    except KeyError:
        # There are no published versions; that's fine!
        latest_id = None

    # Get the ID of the previously uploaded file (if it exists),
    # then delete it so we can upload a new version.
    logger.info("Deleting old file(s)...")
    r = check_status(
        requests.get(
            f"https://{zenodo_url}/api/deposit/depositions/{draft_id}/files",
            params={"access_token": access_token},
        )
    )
    for file in r.json():
        if file["filename"] == file_name:
            FILE_ID = file["id"]
            r = check_status(
                requests.delete(
                    f"https://{zenodo_url}/api/deposit/depositions/{draft_id}/files/{FILE_ID}",
                    params={"access_token": access_token},
                )
            )
            break

    # Upload the new version of the file
    logger.info("Uploading the file...")
    try:
        subprocess.check_output(
            [
                "curl",
                "--progress-bar",
                "--upload-file",
                os.path.join(file_path, file_name),
                "--request",
                "PUT",
                f"{bucket_url}/{file_name}?access_token={access_token}",
            ]
        )
    except Exception as e:
        # Remove traceback info to hide the access token from the output
        with exceptions.no_traceback():
            msg = str(e).replace(access_token, "*" * len(access_token))
            raise exceptions.ZenodoUploadError(msg)

    # Add some metadata
    logger.info("Adding metadata...")
    if shell_cmd == "":
        description = f"{deposit_description}<br/><br/>Created using <a href='https://github.com/rodluger/showyourwork'>showyourwork</a> from <a href='{repo_url}'>this GitHub repo</a>.<br/>"
    else:
        description = f"{deposit_description}<br/><br/>Created using <a href='https://github.com/rodluger/showyourwork'>showyourwork</a> from <a href='{repo_url}'>this GitHub repo</a> using the following command:<br/><br/><pre><code class='language-bash'>{shell_cmd}</code></pre><br/>"
    data = {
        "metadata": {
            "title": deposit_title,
            "upload_type": "dataset",
            "description": description,
            "creators": [{"name": name} for name in deposit_creators],
        }
    }
    r = check_status(
        requests.put(
            f"https://{zenodo_url}/api/deposit/depositions/{draft_id}",
            params={"access_token": access_token},
            data=json.dumps(data),
            headers={"Content-Type": "application/json"},
        )
    )

    if os.getenv("SHOWYOURWORK_DISABLE_UPLOAD", "false") == "true":
        # Don't publish deposits if we're in development mode
        logger.warning("NOTE: Not publishing deposit (SHOWYOURWORK_DISABLE_UPLOAD).")
        pass
    else:
        # Publish the deposit
        logger.info("Publishing the deposit...")
        try:
            r = check_status(
                requests.post(
                    f"https://{zenodo_url}/api/deposit/depositions/{draft_id}/actions/publish",
                    params={"access_token": access_token},
                )
            )
        except exceptions.ZenodoError as e:

            if (
                "New version's files must differ from all previous versions"
                in e.message
            ):
                logger.info("No change in the deposit's files. Aborting.")
                # Revert to the previous deposit (the latest version) ID
                draft_id = latest_id
            else:
                raise e

    # Store the deposit URL; this points to the VERSION id of the record
    deposit_url = f"https://{zenodo_url}/record/{draft_id}"
    with open(f"{file_path}/{file_name}.zenodo", "w") as f:
        print(deposit_url, file=f)


def remove_bad_dotzenodo_files(files):
    """
    If an upload failed b/c of missing authorization, let's
    delete the `.zenodo` file we use for bookkeeping to
    force a reupload attempt the next time we build the
    article.

    """
    for file in files:
        if Path(file).exists():
            with open(file, "r") as f:
                contents = f.readline()
            if "UNAUTHORIZED" in contents:
                Path(file).unlink()