"""
Main Zenodo interface.

"""
from . import exceptions, paths, git
from .subproc import parse_request
from .logging import get_logger
import requests
import tarfile
import os
import json
from pathlib import Path
import subprocess
import shutil

try:
    import snakemake
except ModuleNotFoundError:
    snakemake = None


# Supported tarball extensions
zip_exts = ["tar", "tar.gz", "zip"]


def require_access_token(method):
    def wrapper(self, *args, **kwargs):
        if self.access_token:
            return method(self, *args, **kwargs)
        else:
            raise exceptions.MissingZenodoAccessToken(self.token_name)

    return wrapper


services = {
    "zenodo": {
        "url": "zenodo.org",
        "doi_prefix": "10.5281/zenodo.",
        "token_name": "ZENODO_TOKEN",
        "path": lambda: paths.user().zenodo,
        "name": "Zenodo",
    },
    "sandbox": {
        "url": "sandbox.zenodo.org",
        "doi_prefix": "10.5072/zenodo.",
        "token_name": "SANDBOX_TOKEN",
        "path": lambda: paths.user().sandbox,
        "name": "Zenodo Sandbox",
    },
}


class Zenodo:
    def __init__(self, doi_or_service, **kwargs):
        # Parse input
        if str(doi_or_service).lower() in services.keys():

            # Create a new draft on the given service
            service = services[doi_or_service]
            self.doi_prefix = service["doi_prefix"]
            self.url = service["url"]
            self.token_name = service["token_name"]
            self.access_token = self._get_access_token()
            self.path = service["path"]
            self.service = service["name"]
            self.doi = self._create(**kwargs)
            self.deposit_id = self.doi.split(self.doi_prefix)[1]

        else:

            # Parse the DOI
            self.doi = doi_or_service
            try:
                for service in services.values():
                    if self.doi.startswith(service["doi_prefix"]):
                        self.doi_prefix = service["doi_prefix"]
                        self.deposit_id = self.doi.split(self.doi_prefix)[1]
                        self.url = service["url"]
                        self.token_name = service["token_name"]
                        self.access_token = self._get_access_token()
                        self.path = service["path"]
                        self.service = service["name"]
                        break
                else:
                    raise Exception
            except Exception as e:
                raise exceptions.InvalidZenodoDOI(self.doi)

    def _get_access_token(self, error_if_missing=False):
        """
        Return the API access token (stored in an environment variable).

        """
        return os.getenv(self.token_name, None)

    def get_id_type(self):
        """
        Returns the type of a Zenodo `id` ("version" or "concept").

        Caches the result locally.

        """
        cache_file = self.path() / f"{self.deposit_id}" / "id_type.txt"

        if cache_file.exists():

            # Restore from cache
            with open(cache_file, "r") as f:
                id_type = f.readline().replace("\n", "")

        else:

            # Try to find a published record (no authentication needed)
            try:
                r = requests.get(
                    f"https://{self.url}/api/records/{self.deposit_id}"
                )
                data = r.json()
            except Exception as e:
                r = None
                data = {"status": "", "message": str(e)}

            if (r is None) or (r.status_code > 204):

                # Either is private or doesn't exist.
                # In any event, don't cache it, as this could change.
                return "unknown"

            else:

                # This is a public record
                if int(self.deposit_id) == int(data["conceptrecid"]):
                    id_type = "concept"
                elif int(self.deposit_id) == int(data["id"]):
                    id_type = "version"
                else:
                    id_type = "unknown"

            # Cache it
            cache_file.parents[0].mkdir(exist_ok=True)
            with open(cache_file, "w") as f:
                print(id_type, file=f)

        return id_type

    def _create(self, slug=None, branch=None):
        """
        Create a draft of a Zenodo deposit for the current repo & branch.

        """
        logger = get_logger()
        logger.info(f"Creating a draft deposit on {self.service}...")

        # Get the deposit title
        if branch is None:
            branch = git.get_repo_branch()
        if slug is None:
            slug = git.get_repo_slug()
        title = f"Data for {slug} [{branch}]"

        # Create the draft
        data = parse_request(
            requests.post(
                f"https://{self.url}/api/deposit/depositions",
                params={
                    "access_token": self.access_token,
                },
                json={},
            )
        )

        # Add some minimal metadata
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
        data = parse_request(
            requests.put(
                data["links"]["latest_draft"],
                params={"access_token": self.access_token},
                data=json.dumps(metadata),
                headers={"Content-Type": "application/json"},
            )
        )
        doi = f"{self.doi_prefix}{data['conceptrecid']}"
        logger.info(f"Created draft with concept DOI {doi} on {self.service}.")
        return doi

    @require_access_token
    def upload_file_to_draft(self, draft, file, rule_name, tarball=False):
        """
        Upload a file to a Zenodo draft. Delete the current file
        produced by the same rule, if present.

        """
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
            data = parse_request(
                requests.get(
                    files_url,
                    params={"access_token": self.access_token},
                )
            )
            for entry in data:
                if entry["filename"] == rule_name:
                    file_id = entry["id"]
                    parse_request(
                        requests.delete(
                            f"{files_url}/{file_id}",
                            params={"access_token": self.access_token},
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
        try:
            res = subprocess.run(
                [
                    "curl",
                    "-f",
                    *progress_bar,
                    "-o",
                    "/dev/null",
                    "--upload-file",
                    str(file_to_upload),
                    "--request",
                    "PUT",
                    f"{bucket_url}/{rule_name}?access_token={self.access_token}",
                ]
            )
        except:
            raise exceptions.ZenodoUploadError()

        # Delete the tarball if we created it
        if tarball:
            file_to_upload.unlink()

        # Update the provenance
        rule_hashes[rule_name] = file.name
        metadata["notes"] = json.dumps(rule_hashes, indent=4)
        parse_request(
            requests.put(
                draft["links"]["latest_draft"],
                params={"access_token": self.access_token},
                data=json.dumps({"metadata": metadata}),
                headers={"Content-Type": "application/json"},
            )
        )

    @require_access_token
    def download_file_from_draft(self, draft, file, rule_name, tarball=False):
        """
        Downloads a file from a Zenodo draft.

        """
        # Logger
        logger = get_logger()

        # Get rule hashes for the files currently on Zenodo
        metadata = draft["metadata"]
        notes = metadata.get("notes", "{}")
        try:
            rule_hashes = json.loads(notes)
        except json.JSONDecodeError:
            raise exceptions.InvalidZenodoNotesField()

        # Get the files currently on the remote
        data = parse_request(
            requests.get(
                draft["links"]["files"],
                params={"access_token": self.access_token},
            )
        )

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
                try:
                    res = subprocess.run(
                        [
                            "curl",
                            "-f",
                            f"{url}?access_token={self.access_token}",
                            *progress_bar,
                            "--output",
                            str(file),
                        ]
                    )
                except:
                    raise exceptions.ZenodoDownloadError()

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

    def download_file_from_record(
        self, record, file, rule_name, tarball=False
    ):
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
                subprocess.run(
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

    @require_access_token
    def delete(self):
        """
        Deletes the draft associated with the given concept ID on Zenodo.

        """
        # Logger
        logger = get_logger()

        # Grab the version id
        logger.info(
            f"Deleting {self.service} deposit with concept DOI {self.doi}..."
        )
        r = requests.get(
            f"https://{self.url}/api/deposit/depositions",
            params={
                "q": f"conceptrecid:{self.deposit_id}",
                "all_versions": 1,
                "access_token": self.access_token,
            },
        )
        try:
            for data in r.json():
                if not data["submitted"]:
                    break
            else:
                raise Exception
        except:
            raise exceptions.ZenodoRecordNotFound(self.deposit_id)
        version_id = data["id"]
        parse_request(
            requests.delete(
                f"https://{self.url}/api/deposit/depositions/{version_id}",
                params={
                    "access_token": self.access_token,
                },
            )
        )
        logger.info(f"Successfully deleted deposit {self.doi}.")

    @require_access_token
    def publish(self):
        """
        Publishes the draft associated with the given concept ID on Zenodo.

        """
        # Logger
        logger = get_logger()

        # Grab the version id
        logger.info(
            f"Publishing {self.service} deposit with concept DOI {self.doi}..."
        )
        r = requests.get(
            f"https://{self.url}/api/deposit/depositions",
            params={
                "q": f"conceptrecid:{self.deposit_id}",
                "all_versions": 1,
                "access_token": self.access_token,
            },
        )
        try:
            for data in r.json():
                if not data["submitted"]:
                    break
            else:
                raise Exception
        except:
            raise exceptions.ZenodoRecordNotFound(self.deposit_id)
        version_id = data["id"]
        parse_request(
            requests.post(
                f"https://{self.url}/api/deposit/depositions/{version_id}/actions/publish",
                params={
                    "access_token": self.access_token,
                },
            )
        )
        logger.info(f"Successfully published deposit {self.doi}.")

    def download_file(self, file, rule_name, tarball=False):
        # Logger
        logger = get_logger()

        # Check if there's a draft (and the user has access), and check for
        # a file match. If not, check for existing published versions, and
        # check for a match in each one. If no file is found, raise a cache
        # miss exception, which is caught in the enclosing scope.

        # Check for an existing draft; if found, check for the file in
        # that draft and return if found.
        concept_id = self.deposit_id
        logger.debug(
            f"Attempting to access {self.service} deposit with DOI {self.doi}..."
        )
        r = requests.get(
            f"https://{self.url}/api/deposit/depositions",
            params={
                "q": f"conceptrecid:{concept_id}",
                "all_versions": 1,
                "access_token": self.access_token,
            },
        )
        if r.status_code <= 204:
            try:
                data = r.json()
            except:
                data = []
            if len(data):
                data = data[0]
                draft_url = data.get("links", {}).get("latest_draft", None)
                if draft_url:
                    r = requests.get(
                        draft_url,
                        params={"access_token": self.access_token},
                    )
                    if r.status_code <= 204:
                        try:
                            draft = r.json()
                            self.download_file_from_draft(
                                draft, file, rule_name, tarball=tarball
                            )
                        except exceptions.FileNotFoundOnZenodo:
                            exceptions.restore_trace()
                            logger.debug(
                                f"File {rule_name} not found in deposit with DOI {self.doi}."
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
                logger.debug(
                    f"Failed to access {self.service} deposit with DOI {self.doi}."
                )
        else:
            logger.debug(
                f"Failed to access {self.service} deposit with DOI {self.doi}."
            )
            try:
                data = r.json()
            except:
                pass
            else:
                logger.debug(data["message"])

        # Check for a published record
        logger.debug(
            f"Attempting to access {self.service} record with DOI {self.doi}..."
        )
        r = requests.get(f"https://{self.url}/api/records/{concept_id}")
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
                        "message",
                        f"An error occurred while accessing {self.service}.",
                    ),
                )
        else:
            # There's a published record. Let's search all versions for
            # a file match.
            r = requests.get(
                f"https://{self.url}/api/records",
                params={
                    "q": f'conceptdoi:"{self.doi_prefix}{concept_id}"',
                    "access_token": self.access_token,
                    "all_versions": 1,
                },
            )
            if r.status_code <= 204:
                try:
                    records = r.json().get("hits", {}).get("hits", [])
                except:
                    records = []
                    logger.debug(
                        f"File {rule_name} not found in record with DOI {self.doi}."
                    )
                for record in records[::-1]:
                    try:
                        self.download_file_from_record(
                            record, file, rule_name, tarball=tarball
                        )
                    except exceptions.FileNotFoundOnZenodo:
                        exceptions.restore_trace()
                        logger.debug(
                            f"File {rule_name} not found in record with DOI {self.doi}."
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
                        "message",
                        f"An error occurred while accessing {self.service}.",
                    ),
                )

        # This is caught in the enclosing scope and treated as a cache miss
        raise exceptions.FileNotFoundOnZenodo(file.name)

    @require_access_token
    def upload_file(self, file, rule_name, tarball=False):
        # Logger
        logger = get_logger()

        # Check if a draft already exists, and create it if not.
        # If authentication fails, return with a gentle warning
        concept_id = self.deposit_id
        r = requests.get(
            f"https://{self.url}/api/deposit/depositions",
            params={
                "q": f"conceptrecid:{concept_id}",
                "all_versions": 1,
                "access_token": self.access_token,
            },
        )
        if r.status_code > 204:
            logger.warn(
                f"{self.service} authentication failed. Unable to upload cache for rule {rule_name}."
            )
            try:
                data = r.json()
            except:
                pass
            else:
                logger.debug(data["message"])
            return

        try:
            data = r.json()
        except:
            data = []
        if len(data):
            data = data[0]
        else:
            logger.warn(
                f"{self.service} authentication failed. Unable to upload cache for rule {rule_name}."
            )
            return
        draft_url = data.get("links", {}).get("latest_draft", None)
        if not draft_url and not data["submitted"]:
            draft_url = data["links"]["self"]
        if draft_url:

            # Draft exists
            draft = parse_request(
                requests.get(
                    draft_url,
                    params={"access_token": self.access_token},
                )
            )

        else:

            # Create a new draft
            data = parse_request(
                requests.post(
                    f"https://{self.url}/api/deposit/depositions/{data['id']}/actions/newversion",
                    params={"access_token": self.access_token},
                )
            )
            draft_url = data["links"]["latest_draft"]
            draft = parse_request(
                requests.get(
                    draft_url,
                    params={"access_token": self.access_token},
                )
            )

        self.upload_file_to_draft(draft, file, rule_name, tarball=tarball)

    @require_access_token
    def _download_latest_draft(self):
        # Logger
        logger = get_logger()

        # Grab the deposit
        concept_id = self.deposit_id
        logger.debug(
            f"Attempting to access {self.service} deposit with DOI {self.doi}..."
        )
        r = requests.get(
            f"https://{self.url}/api/deposit/depositions",
            params={
                "q": f"conceptrecid:{concept_id}",
                "all_versions": 1,
                "access_token": self.access_token,
            },
        )
        if r.status_code <= 204:
            try:
                data = r.json()
            except:
                raise exceptions.ZenodoError(
                    message=f"Error accessing latest draft for DOI {self.doi}."
                )

            # Look for a draft
            if len(data):
                data = data[0]
                draft_url = data.get("links", {}).get("latest_draft", None)
                if not draft_url and not data["submitted"]:
                    draft_url = data["links"]["self"]

                # Create a new draft if needed
                if not draft_url:
                    r = requests.post(
                        f"https://{self.url}/api/deposit/depositions/{data['id']}/actions/newversion",
                        params={"access_token": self.access_token},
                    )
                    try:
                        data = r.json()
                    except:
                        raise exceptions.ZenodoError(
                            message=f"Error accessing latest draft for DOI {self.doi}."
                        )
                    draft_url = data["links"]["latest_draft"]

                # Grab the draft
                r = requests.get(
                    draft_url,
                    params={"access_token": self.access_token},
                )
                if r.status_code <= 204:
                    draft = r.json()
                else:
                    raise exceptions.ZenodoError(
                        message=f"Error accessing latest draft for DOI {self.doi}."
                    )
            else:
                raise exceptions.ZenodoError(
                    message=f"Error accessing latest draft for DOI {self.doi}."
                )
        else:
            raise exceptions.ZenodoError(
                message=f"Error accessing latest draft for DOI {self.doi}."
            )

        # Local folder to save to
        cache_folder = self.path() / f"{self.deposit_id}" / "download"
        if cache_folder.exists():
            shutil.rmtree(cache_folder)
        cache_folder.mkdir(exist_ok=True, parents=True)

        # Get metadata & store on disk
        metadata = draft["metadata"]
        with open(cache_folder / ".metadata.json", "w") as f:
            json.dump(metadata, f)

        # Download all files
        data = parse_request(
            requests.get(
                draft["links"]["files"],
                params={"access_token": self.access_token},
            )
        )
        for entry in data:
            url = entry["links"]["download"]
            try:
                res = subprocess.run(
                    [
                        "curl",
                        "-f",
                        f"{url}?access_token={self.access_token}",
                        "--progress-bar",
                        "--output",
                        entry["filename"],
                    ],
                    cwd=cache_folder,
                )
            except:
                raise exceptions.ZenodoDownloadError()

        # Return path to cache folder
        return cache_folder

    @require_access_token
    def copy_draft(self, target_doi_or_service, **kwargs):
        # Logger
        logger = get_logger()
        logger.info(f"Downloading files from {self.doi}...")

        # Download the current draft
        cache_folder = self._download_latest_draft()

        # The target deposit (creates if needed)
        target_deposit = Zenodo(target_doi_or_service, **kwargs)

        # Grab the target deposit
        r = requests.get(
            f"https://{target_deposit.url}/api/deposit/depositions",
            params={
                "q": f"conceptrecid:{target_deposit.deposit_id}",
                "all_versions": 1,
                "access_token": target_deposit.access_token,
            },
        )
        if r.status_code <= 204:
            try:
                data = r.json()
            except:
                raise exceptions.ZenodoError(
                    message=f"Error accessing latest draft for DOI {target_deposit.doi}."
                )

            # Look for a draft
            if len(data):
                data = data[0]
                draft_url = data.get("links", {}).get("latest_draft", None)
                if not draft_url and not data["submitted"]:
                    draft_url = data["links"]["self"]

                # Create a new draft if needed
                if not draft_url:
                    r = requests.post(
                        f"https://{target_deposit.url}/api/deposit/depositions/{data['id']}/actions/newversion",
                        params={"access_token": target_deposit.access_token},
                    )
                    try:
                        data = r.json()
                    except:
                        raise exceptions.ZenodoError(
                            message=f"Error accessing latest draft for DOI {target_deposit.doi}."
                        )
                    draft_url = data["links"]["latest_draft"]

                # Grab the draft
                r = requests.get(
                    draft_url,
                    params={"access_token": target_deposit.access_token},
                )
                if r.status_code <= 204:
                    draft = r.json()
                else:
                    raise exceptions.ZenodoError(
                        message=f"Error accessing latest draft for DOI {target_deposit.doi}."
                    )
            else:
                raise exceptions.ZenodoError(
                    message=f"Error accessing latest draft for DOI {target_deposit.doi}."
                )
        else:
            raise exceptions.ZenodoError(
                message=f"Error accessing latest draft for DOI {target_deposit.doi}."
            )

        # Get the bucket to upload files to
        bucket_url = draft["links"]["bucket"]

        # Upload metadata
        with open(cache_folder / ".metadata.json", "r") as f:
            metadata = json.load(f)
        metadata = {
            "metadata": {
                "title": metadata["title"],
                "upload_type": "dataset",
                "description": metadata["description"],
                "creators": [{"name": "showyourwork"}],
                "notes": metadata.get("notes", "{}"),
            }
        }
        parse_request(
            requests.put(
                draft["links"]["latest_draft"],
                params={"access_token": target_deposit.access_token},
                data=json.dumps(metadata),
                headers={"Content-Type": "application/json"},
            )
        )

        # Upload files
        logger.info(f"Uploading files to {target_deposit.doi}...")
        for file in cache_folder.glob("*"):

            if file.name == ".metadata.json":
                continue

            try:
                res = subprocess.run(
                    [
                        "curl",
                        "-f",
                        "--progress-bar",
                        "-o",
                        "/dev/null",
                        "--upload-file",
                        file.name,
                        "--request",
                        "PUT",
                        f"{bucket_url}/{file.name}?access_token={target_deposit.access_token}",
                    ],
                    cwd=cache_folder,
                )
            except:
                raise exceptions.ZenodoUploadError()

        # We're done
        logger.info(f"Successfully copied {self.doi} to {target_deposit.doi}.")
        return target_deposit.doi
