"""
Main Zenodo interface.

"""
from .exceptions import ShowyourworkException
from pathlib import Path
import requests
import os
import json
import socket
import subprocess
import sys
from contextlib import contextmanager
import requests.packages.urllib3.util.connection as urllib3_cn


# Force IPv4 for Zenodo connections. This is a temporary hack
# that sometimes greatly speeds up access to the API;
# TODO: we need investigate what's actually happening here!
# See https://stackoverflow.com/a/62599037
# and https://stackoverflow.com/a/46972341
urllib3_cn.allowed_gai_family = lambda: socket.AF_INET


@contextmanager
def no_traceback():
    """
    Sets a custom exception handler for the scope of a 'with' block.

    https://stackoverflow.com/a/40347369
    """
    sys.excepthook = no_traceback_excepthook
    yield
    sys.excepthook = sys.__excepthook__


def no_traceback_excepthook(type, value, traceback):
    """
    Print an exception message without the traceback.

    """
    print(": ".join([str(type.__name__), str(value)]))


def check_status(r):
    """
    Parse a requests return object and raise a custom exception
    for a >200-level status code.

    """
    if r.status_code > 204:
        data = r.json()
        for error in data.get("errors", []):
            data["message"] += " " + error["message"]
        raise ShowyourworkException(
            "Zenodo error {}: {}".format(data["status"], data["message"])
        )
    return r


def get_access_token(token_name, error_if_missing=False):
    """
    Return the access token stored in the environment variable `token_name`.

    """
    access_token = os.getenv(token_name, None)
    if error_if_missing and access_token is None:
        raise ValueError(
            f"Zenodo access token `{token_name}` not found. "
            "This should be set as an environment variable "
            "and/or a repository secret."
        )
    return access_token


def reserve():
    """
    Pre-reserve a concept DOI on Zenodo.
    This is a user-facing command line utility called from the Makefile.
    Not to be run on GitHub Actions!

    """
    # Zenodo Sandbox (for testing) or Zenodo?
    while True:
        try:
            sandbox = input("Pre-reserve on Zenodo Sandbox (y/n) [n]: ")
        except (EOFError, KeyboardInterrupt):
            print("[Interrupted.]")
            return
        if sandbox == "" or sandbox == "n":
            zenodo_url = "zenodo.org"
            break
        elif sandbox == "y":
            zenodo_url = "sandbox.zenodo.org"
            break
        else:
            print("Please enter either 'y' or 'n'.")

    # Default token name
    if sandbox == "y":
        default_token_name = "ZENODO_SANDBOX_TOKEN"
    else:
        default_token_name = "ZENODO_TOKEN"

    # Zenodo access token
    while True:
        try:
            token_name = input(
                f"Name of environment variable containing Zenodo API token [{default_token_name}]: "
            )
        except (EOFError, KeyboardInterrupt):
            print("[Interrupted.]")
            return
        if token_name == "":
            token_name = default_token_name
        access_token = get_access_token(token_name, error_if_missing=False)
        if access_token is None:
            print(f"Access token {token_name} not found.")
        else:
            break

    # Create a new deposit
    print("Pre-reserving Zenodo concept DOI...")
    r = check_status(
        requests.post(
            f"https://{zenodo_url}/api/deposit/depositions",
            params={
                "access_token": access_token,
            },
            json={},
        )
    )

    if r is not None:

        # Get the deposit id
        DEPOSIT_ID = r.json()["conceptrecid"]

        # We're done
        print(f"Deposit reserved on {zenodo_url} with concept DOI id {DEPOSIT_ID}.")


def get_id_type(
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
                    raise ShowyourworkException(
                        f"Record/deposit with id {deposit_id} not found.",
                        context=f"The provided `id` {deposit_id} does "
                        "not seem to be a valid concept or version Zenodo id.",
                        brief=f"Invalid Zenodo id {deposit_id}.",
                    )

        else:

            # Something unexpected happened...
            status = data.get("status", "unknown")
            message = data.get("message", "An error occurred while accessing Zenodo.")
            context = "Zenodo error {}: {}".format(status, message)
            raise ShowyourworkException(
                message,
                context=context,
                brief="An error occurred while accessing Zenodo.",
            )

    else:

        # This is a public record
        if int(deposit_id) == int(data["conceptrecid"]):
            return "concept"
        elif int(deposit_id) == int(data["id"]):
            return "version"
        else:
            raise ShowyourworkException(
                f"Record/deposit with id {deposit_id} not found.",
                context=f"The provided `id` {deposit_id} does "
                "not seem to be a valid concept or version Zenodo id.",
                brief=f"Invalid Zenodo id {deposit_id}.",
            )


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
                    raise ShowyourworkException(
                        f"Record/deposit with concept id {concept_id} not found.",
                        context=f"The provided `id` {concept_id} does "
                        "not seem to be a valid Zenodo concept id.",
                        brief=f"Invalid Zenodo concept id {concept_id}.",
                    )

        except StopIteration:

            # We're good!
            pass

        else:

            # We reached the end and found nothing
            raise ShowyourworkException(
                f"Record/deposit with concept id {concept_id} not found.",
                context=f"The provided `id` {concept_id} does "
                "not seem to be a valid Zenodo concept id.",
                brief=f"Invalid Zenodo concept id {concept_id}.",
            )

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

        raise ShowyourworkException(
            f"Record/deposit with concept id {concept_id} not found.",
            context=f"The provided `id` {concept_id} does "
            "not seem to be a valid Zenodo concept id.",
            brief=f"Invalid Zenodo concept id {concept_id}.",
        )

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
    print("Deleting old file(s)...")
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
    print("Uploading the file...")
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
        msg = str(e).replace(access_token, "*" * len(access_token))
        with no_traceback():
            # Don't display the traceback, which will usually show
            # the command we invoked containing the access token.
            # Hide the access token from the error message.
            raise Exception(msg)
    """
    # Old method using requests.
    # This is safer (no chance of the access token getting printed
    # to the logs) but AFAIK there's no way to add a progress bar.
    # There are lots of solutions online that use `requests_toolbelt`,
    # but those all require us to use the POST (not PUT) method.
    with open(os.path.join(file_path, file_name), "rb") as fp:
        r = check_status(
            requests.put(
                "%s/%s" % (bucket_url, file_name),
                data=fp,
                params={"access_token": access_token},
            )
        )
    """

    # Add some metadata
    print("Adding metadata...")
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
        print("NOTE: Not publishing deposit (SHOWYOURWORK_DISABLE_UPLOAD).")
        pass
    else:
        # Publish the deposit
        print("Publishing the deposit...")
        try:
            r = check_status(
                requests.post(
                    f"https://{zenodo_url}/api/deposit/depositions/{draft_id}/actions/publish",
                    params={"access_token": access_token},
                )
            )
        except ShowyourworkException as e:

            if (
                "New version's files must differ from all previous versions"
                in e.message
            ):
                print("No change in the deposit's files. Aborting.")
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