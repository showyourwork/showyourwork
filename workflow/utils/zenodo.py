"""
Main Zenodo interface.

"""
from . import exceptions, paths, git
from .logging import get_logger
import requests
import tarfile
import os
import json
from pathlib import Path
import snakemake
import subprocess


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
        raise exceptions.ZenodoError(status=data["status"], message=data["message"])
    return r


def get_access_token(token_name="ZENODO_TOKEN", error_if_missing=False):
    """
    Return the access token stored in the environment variable `token_name`.

    """
    access_token = os.getenv(token_name, None)
    if error_if_missing and access_token is None:
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
            raise exceptions.ZenodoError(status=data["status"], message=data["message"])

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
        tmp = paths.zenodo_sandbox
    else:
        tmp = paths.zenodo
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
    for s in ["https://", "http://"]:
        url = url.replace(s, "")
    try:
        branch = snakemake.workflow.config["branch"]
    except:
        branch = git.get_repo_branch()
    return f"Data for {url}@{branch}"


def create_draft():
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
            "https://zenodo.org/api/deposit/depositions",
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
        "Data automatically uploaded by the showyourwork workflow. "
        "Each of the files in this deposit were generated from a user-defined "
        "Snakemake rule with the same name and hash specified in the Notes "
        "section below. "
        "Please visit github.com/rodluger/showyourwork and/or the article source "
        "repository for more information."
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

    return data


def get_draft():
    """
    Get the Zenodo deposit draft for the current repo & branch.
    Creates it if it doesn't exist.

    """
    # Get the Zenodo token
    access_token = get_access_token(error_if_missing=True)

    # Search
    title = get_deposit_title()
    r = requests.get(
        "https://zenodo.org/api/deposit/depositions",
        params={
            "q": f'title:"{title}"',
            "access_token": access_token,
            "status": "draft",
            "size": 10,
            "page": 1,
            "sort": "mostrecent",
            "all_versions": 0,
        },
    )
    data = r.json()
    if len(data):

        # We found a matching draft; return the most recent one.
        r = check_status(
            requests.get(
                data[0]["links"]["latest_draft"],
                params={"access_token": access_token},
            )
        )
        return r.json()

    else:

        # We'll create a new draft
        return create_draft()


def upload_file_to_draft(file, rule_name, tarball=False):
    """
    Upload a file to a Zenodo draft. Delete the current file
    produced by the same rule, if present.

    Note that the name of the file _is_ the hash of the rule
    that generated it (courtesy of Snakemake).

    """
    # Get the Zenodo token
    access_token = get_access_token(error_if_missing=True)

    # Get the current draft
    draft = get_draft()

    # Get rule hashes for the files currently on Zenodo
    metadata = draft["metadata"]
    notes = metadata.get("notes", "{}")
    try:
        rule_hashes = json.loads(notes)
    except json.JSONDecodeError:
        # TODO
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
    try:
        progress_bar = (
            ["--progress-bar"]
            if not snakemake.workflow.config["github_actions"]
            else []
        )
        subprocess.check_output(
            [
                "curl",
                *progress_bar,
                "--upload-file",
                file_to_upload,
                "--request",
                "PUT",
                f"{bucket_url}/{rule_name}?access_token={access_token}",
            ]
        )
    except Exception as e:
        msg = str(e).replace(access_token, "*" * len(access_token))
        with exceptions.no_traceback():
            # Don't display the traceback, which will usually show
            # the command we invoked containing the access token.
            # Hide the access token from the error message.
            raise Exception(msg)

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


def download_file_from_draft(file, rule_name, tarball=False):
    """
    Downloads a file from a Zenodo draft.

    """
    # Get the Zenodo token
    access_token = get_access_token(error_if_missing=True)

    # Get the current draft
    draft = get_draft()

    # Get rule hashes for the files currently on Zenodo
    metadata = draft["metadata"]
    notes = metadata.get("notes", "{}")
    try:
        rule_hashes = json.loads(notes)
    except json.JSONDecodeError:
        # TODO
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
            try:
                progress_bar = (
                    ["--progress-bar"]
                    if not snakemake.workflow.config["github_actions"]
                    else []
                )
                subprocess.check_output(
                    [
                        "curl",
                        f"{url}?access_token={access_token}",
                        *progress_bar,
                        "--output",
                        str(file),
                    ]
                )
            except Exception as e:
                msg = str(e).replace(access_token, "*" * len(access_token))
                with exceptions.no_traceback():
                    # Don't display the traceback, which will usually show
                    # the command we invoked containing the access token.
                    # Hide the access token from the error message.
                    raise Exception(msg)

            # If it's a directory tarball, extract it
            if tarball:
                os.rename(file, f"{file}.tar.gz")
                with tarfile.open(f"{file}.tar.gz") as tb:
                    tb.extractall(file)
                Path(f"{file}.tar.gz").unlink()

            return

    # This is caught in the enclosing scope and treated as a cache miss
    raise exceptions.FileNotFoundOnZenodo(file.name)
