"""
Main Zenodo interface.

"""
from . import exceptions, paths, git
from .subproc import run
from .logging import get_logger
import requests
import tarfile
import os
import json
from pathlib import Path
import subprocess

try:
    import snakemake
except ModuleNotFoundError:
    snakemake = None


# Zenodo base URLs
zenodo_url = {"zenodo": "zenodo.org", "zenodo_sandbox": "sandbox.zenodo.org"}


# Supported tarball extensions
zip_exts = ["tar.gz"]


def check_status(r):
    """
    Parse a requests return object and raise a custom exception
    for a >200-level status code.

    """
    if r.status_code > 204:
        try:
            data = r.json()
        except:
            raise exceptions.ZenodoError(status=r.status_code)
        for error in data.get("errors", []):
            data["message"] += " " + error["message"]
        raise exceptions.ZenodoError(
            status=data["status"], message=data["message"]
        )
    return r


def get_access_token(token_name="ZENODO_TOKEN", error_if_missing=False):
    """
    Return the access token stored in the environment variable `token_name`.

    """
    access_token = os.getenv(token_name, None)
    if error_if_missing and not access_token:
        raise exceptions.MissingZenodoAccessToken(token_name)
    return access_token


def _get_id_type(deposit_id, zenodo_url="zenodo.org"):
    """
    Determines whether a given Zenodo `id` corresponds to
    a concept id or a version id.

    """
    # Try to find a published record (no authentication needed)
    r = requests.get(f"https://{zenodo_url}/api/records/{deposit_id}")
    data = r.json()

    if r.status_code > 204:

        if "PID is not registered" in data.get("message", ""):

            # No published records found
            raise exceptions.ZenodoRecordNotFound(deposit_id)

        else:

            # Something unexpected happened...
            raise exceptions.ZenodoError(
                status=data["status"], message=data["message"]
            )

    else:

        # This is a public record
        if int(deposit_id) == int(data["conceptrecid"]):
            return "concept"
        elif int(deposit_id) == int(data["id"]):
            return "version"
        else:
            raise exceptions.ZenodoRecordNotFound(deposit_id)


def get_id_type(deposit_id, zenodo_url="zenodo.org"):
    """
    Returns the type of a Zenodo `id` ("version" or "concept").

    Caches the result locally.

    """
    if "sandbox" in zenodo_url:
        tmp = paths.user().zenodo_sandbox
    else:
        tmp = paths.user().zenodo
    cache_file = tmp / f"{deposit_id}" / "id_type.txt"

    if cache_file.exists():

        with open(cache_file, "r") as f:
            id_type = f.readline().replace("\n", "")

    else:

        cache_file.parents[0].mkdir(exist_ok=True)
        id_type = _get_id_type(deposit_id, zenodo_url=zenodo_url)
        with open(cache_file, "w") as f:
            print(id_type, file=f)

    return id_type


def get_deposit_title():
    """
    Return the title of the Zenodo deposit for the current repo branch.

    """
    # Repo info
    try:
        url = snakemake.workflow.config["url"]
    except:
        url = git.get_repo_url()
    for s in [
        "https://github.com/",
        "http://github.com/",
        "https://",
        "http://",
    ]:
        url = url.replace(s, "")
    try:
        branch = snakemake.workflow.config["branch"]
    except:
        branch = git.get_repo_branch()
    return f"Data for {url} [{branch}]"


def create_deposit():
    """
    Create a draft of a Zenodo deposit for the current repo & branch.

    """
    logger = get_logger()
    logger.info("Creating a draft deposit on Zenodo...")

    # Get the Zenodo token
    access_token = get_access_token(error_if_missing=True)

    # Create the draft
    r = check_status(
        requests.post(
            f"https://zenodo.org/api/deposit/depositions",
            params={
                "access_token": access_token,
            },
            json={},
        )
    )

    # Add some minimal metadata
    data = r.json()
    title = get_deposit_title()
    description = (
        "Data automatically uploaded by the <code>showyourwork</code> workflow. "
        "Each of the files in this deposit were generated from a user-defined "
        "<code>Snakemake</code> rule with the same name and hash specified "
        "below. Please visit "
        "<a href='https://github.com/showyourwork'>github.com/showyourwork</a> "
        "and/or the article source repository for more information."
    )
    metadata = {
        "metadata": {
            "title": title,
            "upload_type": "dataset",
            "description": description,
            "creators": [{"name": "showyourwork"}],
            "notes": "{}",
        }
    }
    r = check_status(
        requests.put(
            data["links"]["latest_draft"],
            params={"access_token": access_token},
            data=json.dumps(metadata),
            headers={"Content-Type": "application/json"},
        )
    )
    data = r.json()
    logger.info(f"Created draft with id {data['id']} on Zenodo.")

    # Return the ids
    concept_id = data["conceptrecid"]
    return concept_id


def upload_file_to_draft(draft, file, rule_name, tarball=False):
    """
    Upload a file to a Zenodo draft. Delete the current file
    produced by the same rule, if present.

    """
    # Get the Zenodo token
    access_token = get_access_token(error_if_missing=True)

    # Get rule hashes for the files currently on Zenodo
    metadata = draft["metadata"]
    notes = metadata.get("notes", "{}")
    try:
        rule_hashes = json.loads(notes)
    except json.JSONDecodeError:
        raise exceptions.InvalidZenodoNotesField()

    # Search for an existing file on Zenodo
    rule_hash_on_zenodo = rule_hashes.get(rule_name, None)
    if rule_hash_on_zenodo == file.name:
        # The file is up to date
        return
    elif rule_hash_on_zenodo:
        # Delete the existing file
        files_url = draft["links"]["files"]
        r = check_status(
            requests.get(
                files_url,
                params={"access_token": access_token},
            )
        )
        data = r.json()
        for entry in data:
            if entry["filename"] == rule_name:
                file_id = entry["id"]
                r = check_status(
                    requests.delete(
                        f"{files_url}/{file_id}",
                        params={"access_token": access_token},
                    )
                )
                break

    # If it's a directory, tar it up
    if tarball:
        with tarfile.open(f"{file}.tar.gz", "w:gz") as tb:
            tb.add(file, arcname=".")
        file_to_upload = Path(f"{file}.tar.gz")
    else:
        file_to_upload = file

    # Use curl to upload the file so we have a progress bar
    bucket_url = draft["links"]["bucket"]
    progress_bar = (
        ["--progress-bar"]
        if not snakemake.workflow.config["github_actions"]
        else []
    )
    run(
        [
            "curl",
            *progress_bar,
            "--upload-file",
            file_to_upload,
            "--request",
            "PUT",
            f"{bucket_url}/{rule_name}?access_token={access_token}",
        ],
        secrets=[access_token],
    )

    # Delete the tarball if we created it
    if tarball:
        file_to_upload.unlink()

    # Update the provenance
    rule_hashes[rule_name] = file.name
    metadata["notes"] = json.dumps(rule_hashes, indent=4)
    r = check_status(
        requests.put(
            draft["links"]["latest_draft"],
            params={"access_token": access_token},
            data=json.dumps({"metadata": metadata}),
            headers={"Content-Type": "application/json"},
        )
    )


def download_file_from_draft(draft, file, rule_name, tarball=False):
    """
    Downloads a file from a Zenodo draft.

    """
    # Logger
    logger = get_logger()

    # Get the Zenodo token
    access_token = get_access_token(error_if_missing=False)
    if not access_token:
        logger.debug("Zenodo access token not provided; can't search drafts.")
        exceptions.FileNotFoundOnZenodo(rule_name)

    # Get rule hashes for the files currently on Zenodo
    metadata = draft["metadata"]
    notes = metadata.get("notes", "{}")
    try:
        rule_hashes = json.loads(notes)
    except json.JSONDecodeError:
        raise exceptions.InvalidZenodoNotesField()

    # Get the files currently on the remote
    r = check_status(
        requests.get(
            draft["links"]["files"],
            params={"access_token": access_token},
        )
    )
    data = r.json()

    # Look for a match
    for entry in data:
        if (
            entry["filename"] == rule_name
            and rule_hashes.get(rule_name, None) == file.name
        ):

            # Download it
            url = entry["links"]["download"]
            progress_bar = (
                ["--progress-bar"]
                if not snakemake.workflow.config["github_actions"]
                else []
            )
            run(
                [
                    "curl",
                    f"{url}?access_token={access_token}",
                    *progress_bar,
                    "--output",
                    str(file),
                ],
                secrets=[access_token],
            )

            # If it's a directory tarball, extract it
            if tarball:
                os.rename(file, f"{file}.tar.gz")
                with tarfile.open(f"{file}.tar.gz") as tb:
                    tb.extractall(file)
                Path(f"{file}.tar.gz").unlink()

            return

        elif entry["filename"] == rule_name:

            logger.debug(
                f"File {rule_name} found, but it has the wrong hash. Skipping..."
            )
            break

    # This is caught in the enclosing scope and treated as a cache miss
    raise exceptions.FileNotFoundOnZenodo(rule_name)


def download_file_from_record(record, file, rule_name, tarball=False):
    """
    Downloads a file from a published Zenodo record.

    """
    # Get rule hashes for the files currently on Zenodo
    metadata = record["metadata"]
    notes = metadata.get("notes", "{}")
    try:
        rule_hashes = json.loads(notes)
    except json.JSONDecodeError:
        raise exceptions.InvalidZenodoNotesField()

    # Look for a match
    for entry in record["files"]:
        if (
            entry["key"] == rule_name
            and rule_hashes.get(rule_name, None) == file.name
        ):

            # Download it
            url = entry["links"]["self"]
            progress_bar = (
                ["--progress-bar"]
                if not snakemake.workflow.config["github_actions"]
                else []
            )
            run(
                [
                    "curl",
                    url,
                    *progress_bar,
                    "--output",
                    str(file),
                ]
            )

            # If it's a directory tarball, extract it
            if tarball:
                os.rename(file, f"{file}.tar.gz")
                with tarfile.open(f"{file}.tar.gz") as tb:
                    tb.extractall(file)
                Path(f"{file}.tar.gz").unlink()

            return

        elif entry["key"] == rule_name:

            logger.debug(
                f"File {rule_name} found, but it has the wrong hash. Skipping..."
            )
            break

    # This is caught in the enclosing scope and treated as a cache miss
    raise exceptions.FileNotFoundOnZenodo(rule_name)


def download_file_from_zenodo(file, rule_name, concept_id, tarball=False):
    # Logger
    logger = get_logger()

    # Check if there's a draft (and the user has access), and check for
    # a file match. If not, check for existing published versions, and
    # check for a match in each one. If no file is found, raise a cache
    # miss exception, which is caught in the enclosing scope.

    # Get the Zenodo token
    access_token = get_access_token(error_if_missing=False)

    # Check for an existing draft; if found, check for the file in
    # that draft and return if found.
    logger.debug(f"Attempting to access Zenodo deposit {concept_id}...")
    r = requests.get(
        f"https://zenodo.org/api/deposit/depositions",
        params={
            "q": f"conceptrecid:{concept_id}",
            "all_versions": 1,
            "access_token": access_token,
        },
    )
    if r.status_code <= 204:
        data = r.json()
        if len(data):
            data = data[0]
            draft_url = data.get("links", {}).get("latest_draft", None)
            if draft_url:
                r = requests.get(
                    draft_url,
                    params={"access_token": access_token},
                )
                if r.status_code <= 204:
                    draft = r.json()
                    try:
                        download_file_from_draft(
                            draft, file, rule_name, tarball=tarball
                        )
                    except exceptions.FileNotFoundOnZenodo:
                        logger.debug(
                            f"File {rule_name} not found in deposit {concept_id}."
                        )
                    else:
                        return
                else:
                    logger.debug(
                        f"Something went wrong accessing {draft_url}."
                    )
                    try:
                        data = r.json()
                    except:
                        pass
                    else:
                        logger.debug(data["message"])

    else:
        logger.debug(f"Failed to access Zenodo deposit {concept_id}.")
        try:
            data = r.json()
        except:
            pass
        else:
            logger.debug(data["message"])

    # Check for a published record
    logger.debug(f"Attempting to access Zenodo record {concept_id}...")
    r = requests.get(f"https://zenodo.org/api/records/{concept_id}")
    if r.status_code > 204:
        try:
            data = r.json()
        except:
            data = {}
        if "PID is not registered" in data.get("message", ""):
            # There is no published record with this id
            pass
        else:
            # Something unexpected happened
            raise exceptions.ZenodoError(
                status=data.get("status", "unknown"),
                message=data.get(
                    "message", "An error occurred while accessing Zenodo."
                ),
            )
    else:
        # There's a published record. Let's search all versions for
        # a file match.
        r = requests.get(
            f"https://zenodo.org/api/records",
            params={
                "q": f'conceptdoi:"10.5281/zenodo.{concept_id}"',
                "access_token": access_token,
                "all_versions": 1,
            },
        )
        if r.status_code <= 204:
            records = r.json().get("hits", {}).get("hits", [])
            for record in records[::-1]:
                try:
                    download_file_from_record(
                        record, file, rule_name, tarball=tarball
                    )
                except exceptions.FileNotFoundOnZenodo:
                    logger.debug(
                        f"File {rule_name} not found in record {concept_id}."
                    )
                else:
                    return
        else:
            # Something unexpected happened
            try:
                data = r.json()
            except:
                data = {}
            raise exceptions.ZenodoError(
                status=data.get("status", "unknown"),
                message=data.get(
                    "message", "An error occurred while accessing Zenodo."
                ),
            )

    # This is caught in the enclosing scope and treated as a cache miss
    raise exceptions.FileNotFoundOnZenodo(file.name)


def upload_file_to_zenodo(file, rule_name, concept_id, tarball=False):
    # Logger
    logger = get_logger()

    # Get the Zenodo token
    access_token = get_access_token(error_if_missing=False)
    if not access_token:
        logger.warn(
            f"Zenodo access token not provided. Unable to upload cache for rule {rule_name}."
        )
        return

    # Check if a draft already exists, and create it if not.
    # If authentication fails, return with a gentle warning
    r = requests.get(
        f"https://zenodo.org/api/deposit/depositions",
        params={
            "q": f"conceptrecid:{concept_id}",
            "all_versions": 1,
            "access_token": access_token,
        },
    )
    if r.status_code > 204:
        logger.warn(
            f"Zenodo authentication failed. Unable to upload cache for rule {rule_name}."
        )
        try:
            data = r.json()
        except:
            pass
        else:
            logger.debug(data["message"])
        return

    data = r.json()
    if len(data):
        data = data[0]
    else:
        logger.warn(
            f"Zenodo authentication failed. Unable to upload cache for rule {rule_name}."
        )
        return
    draft_url = data.get("links", {}).get("latest_draft", None)
    if draft_url:

        # Draft exists
        r = check_status(
            requests.get(
                draft_url,
                params={"access_token": access_token},
            )
        )
        draft = r.json()

    else:

        # Create a new draft
        r = check_status(
            requests.post(
                data["links"]["newversion"],
                params={"access_token": access_token},
            )
        )
        draft_url = r.json()["links"]["latest_draft"]
        r = check_status(
            requests.get(
                draft_url,
                params={"access_token": access_token},
            )
        )
        draft = r.json()

    upload_file_to_draft(draft, file, rule_name, tarball=tarball)
