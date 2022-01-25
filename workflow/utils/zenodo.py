"""
Main Zenodo interface.

"""
from . import exceptions, paths
from .config import get_snakemake_variable
import requests
import os
import json
import hashlib
from pathlib import Path
import snakemake


# Default evironment variable names
default_token_name = {
    "zenodo": "ZENODO_TOKEN",
    "zenodo_sandbox": "ZENODO_SANDBOX_TOKEN",
}

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


def get_access_token(token_name, error_if_missing=False):
    """
    Return the access token stored in the environment variable `token_name`.

    """
    access_token = os.getenv(token_name, None)
    if error_if_missing and access_token is None:
        raise exceptions.MissingZenodoAccessToken(token_name)
    return access_token


def _get_user_id(zenodo_url="zenodo.org", token_name="ZENODO_TOKEN"):
    """
    Return the internal user ID associated with a Zenodo API token.

    """
    # Get the access token
    try:
        access_token = get_access_token(token_name, error_if_missing=True)
    except exceptions.MissingZenodoAccessToken:
        return None

    # Create a new deposit so we can grab the user ID
    try:
        r = check_status(
            requests.post(
                f"https://{zenodo_url}/api/deposit/depositions",
                params={"access_token": access_token},
                data="{}",
                headers={"Content-Type": "application/json"},
            )
        )
    except exceptions.ZenodoError:
        return None

    data = r.json()
    deposit_id = data.get("id", None)
    user_id = data.get("owner", None)
    if user_id:
        user_id = int(user_id)

    # Delete the deposit
    if deposit_id:
        try:
            r = check_status(
                requests.delete(
                    f"https://{zenodo_url}/api/deposit/depositions/{deposit_id}",
                    params={"access_token": access_token},
                )
            )
        except exceptions.ZenodoError:
            # It's... probably fine
            pass

    return user_id


def get_user_id(zenodo_url="zenodo.org", token_name="ZENODO_TOKEN"):
    """
    Return the internal user ID associated with a Zenodo API token.

    A few notes:

        - Caches the result locally if successful.
        - We DO NOT store the API token on disk, but instead its
          SHA-256 hash.
        - If the user is not authenticated (missing or bad API token),
          returns `None`.

    """
    if "sandbox" in zenodo_url:
        tmp = paths.zenodo_sandbox_ids
    else:
        tmp = paths.zenodo_ids

    # Get the access token
    try:
        access_token = get_access_token(token_name, error_if_missing=True)
    except exceptions.MissingZenodoAccessToken:
        return None

    # Unique (but secure) cache file for the provided access token
    access_token_hash = hashlib.sha256(access_token.encode()).hexdigest()
    cache_file = tmp / f"{access_token_hash}.txt"

    if cache_file.exists():

        with open(cache_file, "r") as f:
            user_id = int(f.readline().replace("\n", ""))

    else:

        user_id = _get_user_id(zenodo_url=zenodo_url, token_name=token_name)
        if user_id:
            with open(cache_file, "w") as f:
                print(user_id, file=f)

    return user_id


def _get_owner_ids(data):
    """
    Infer the owner IDs given a data dict returned by a GET deposition request.

    """
    if "owner" in data:
        return [int(data["owner"])]
    elif "owners" in data:
        return [int(i) for i in data["owners"]]
    else:
        # TODO
        raise exceptions.ZenodoError()


def _get_id_info(
    deposit_id,
    zenodo_url="zenodo.org",
    token_name="ZENODO_TOKEN",
    search_drafts=True,
    max_pages=100,
    results_per_page=10,
):
    """
    Determines whether a given Zenodo `id` corresponds to
    a concept id or a version id (and raises an error otherwise).

    Returns a tuple: the record id type and the ID(s) of the owner(s).

    """
    # Try to find a published record (no authentication needed)
    r = requests.get(f"https://{zenodo_url}/api/records/{deposit_id}")
    data = r.json()

    if r.status_code > 204:

        if "PID is not registered" in data.get("message", ""):

            # No published records found
            if not search_drafts:
                raise exceptions.ZenodoRecordNotFound(deposit_id)

            # Get the access token
            access_token = get_access_token(token_name, error_if_missing=True)

            # Search for a draft deposit with this concept or
            # version ID (authentication needed)
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
                        return "concept", _get_owner_ids(deposit)

                    elif int(deposit["id"]) == deposit_id:

                        # We found a deposit draft with this version id
                        return "version", _get_owner_ids(deposit)

                if len(data) == 0:

                    # We reached the end and found nothing
                    raise exceptions.ZenodoRecordNotFound(deposit_id)

        else:

            # Something unexpected happened...
            raise exceptions.ZenodoError(status=data["status"], message=data["message"])

    else:

        # This is a public record
        if int(deposit_id) == int(data["conceptrecid"]):
            return "concept", _get_owner_ids(data)
        elif int(deposit_id) == int(data["id"]):
            return "version", _get_owner_ids(data)
        else:
            raise exceptions.ZenodoRecordNotFound(deposit_id)


def get_id_info(deposit_id, zenodo_url="zenodo.org", **kwargs):
    """
    Determines whether a given Zenodo `id` corresponds to
    a concept id or a version id (and raises an error otherwise).

    Returns a tuple: the record id type and a list of owner IDs.

    Caches the result locally.

    """
    if "sandbox" in zenodo_url:
        tmp = paths.zenodo_sandbox
    else:
        tmp = paths.zenodo
    cache_file = tmp / f"{deposit_id}" / "info.json"

    if cache_file.exists():

        with open(cache_file, "r") as f:
            data = json.load(f)
            id_type = data["id_type"]
            owner_ids = data["owner_ids"]

    else:

        cache_file.parents[0].mkdir(exist_ok=True)
        id_type, owner_ids = _get_id_info(deposit_id, zenodo_url=zenodo_url, **kwargs)
        data = {"id_type": id_type, "owner_ids": owner_ids}
        with open(cache_file, "w") as f:
            json.dump(data, f)

    return id_type, owner_ids


class _RuleHash:
    def __init__(self, rule):
        self.output = [Path(file).as_posix() for file in rule.output]

    def _get_direct_dependencies(self, file):
        """
        Get the inputs to all jobs that produce `file`.

        """
        # Note that there should only be at most 1 job if
        # the workflow has no ambiguous rules (and even if it
        # does, the DAG should choose a single job, I think)
        inputs = [job.input for job in self.dag.jobs if file in job.rule.output]

        # Return a set
        return set([Path(item).as_posix() for sublist in inputs for item in sublist])

    def _get_all_dependencies(self, file):
        """
        Recursively get all dependencies of a file.

        """
        deps = self._get_direct_dependencies(file)
        new_deps = set()
        for dep in deps:
            new_deps.update(self._get_all_dependencies(dep))
        deps.update(new_deps)
        return deps

    def _get_file_metadata(self, file, buf_sz=2 ** 16):
        """
        Return metadata for a dependency.

        Dependencies can either be "static" (e.g., files comitted to the repo)
        or they can be programmatically generated by a job. In the former case,
        we return the file's SHA256 hash. In the latter case, we return the
        SHA256 of the generating rule (since the file may not exist yet).

        """
        for job in self.dag.jobs:

            if file in job.rule.output:

                # Get all inputs, including script & conda env
                inputs = set(job.unique_input)
                if job.rule.conda_env:
                    inputs.add(job.rule.conda_env)
                if job.rule.script:
                    inputs.add(job.rule.script)
                inputs = sorted(inputs)

                # Assemble all relevant metadata for the rule
                # Note that the `script`, `notebook`, and `conda` files
                # are automatically included in `inputs` for all user-defined
                # rules; see `userrules.py`
                # TODO: Record these in a temporary JSON file for debugging
                meta = {
                    "input": inputs,
                    "params": dict(job.params),
                    "shell": job.shellcmd,
                }

                meta_str = json.dumps(meta, sort_keys=True)
                sha256 = hashlib.sha256(meta_str.encode())
                return {file: sha256.hexdigest()}

        else:

            if not (paths.user / file).exists():
                # TODO (this shouldn't happen)
                raise exceptions.ShowyourworkException()

            # Compute the file hash
            sha256 = hashlib.sha256()
            with open(paths.user / file, "rb") as f:
                while True:
                    data = f.read(buf_sz)
                    if not data:
                        break
                    sha256.update(data)
            return {file: sha256.hexdigest()}

    def _get_rule_metadata(self):
        """
        Return a dict containing metadata about all the inputs and
        parameters of a Snakemake rule.

        This dict contains the SHA256 hash of all static inputs and
        the settings for the rules generating all dynamic inputs to
        the Snakemake rule.

        """
        # Wait until the DAG has been built
        snakemake.workflow.checkpoints.syw__dag.get()

        # Get the workflow graph
        self.dag = get_snakemake_variable("dag", None)
        if self.dag is None:
            raise exceptions.ShowyourworkException("Unable to obtain the workflow DAG.")

        # TODO: Handle workflows with user-defined checkpoints. Currently
        # there's no way to infer the true dependencies of a file if they
        # originate from a checkpoint without running the workflow twice!
        if len(list(self.dag.checkpoint_jobs)):
            raise exceptions.NotImplementedError(
                "User-defined checkpoints are not currently supported."
            )

        for file in self.output:

            # Get all dependencies of the file
            deps = self._get_all_dependencies(file)

            # Get metadata for each dependency
            meta = {}
            for dep in sorted(deps):
                meta.update(self._get_file_metadata(dep))

        return meta

    def _get_rule_hash(self):
        """
        Return a unique hash for the Snakemake rule.

        Two rules will have the same hash if and only if their static input
        files have the same SHA256 hash and if their programmatically-generated
        input files were created by rules with the same hash.

        """
        meta = self._get_rule_metadata()
        meta_str = json.dumps(meta, sort_keys=True)
        sha256 = hashlib.sha256(meta_str.encode())
        return sha256.hexdigest()


class disable_if_downloadable(_RuleHash):
    """
    When used as an input file function to a rule, dynamically disables
    a custom user rule if the (exact) output of the rule can be downloaded
    from Zenodo.

    """

    def __call__(self, wildcards):

        print(self._get_rule_hash())

        # TODO: Insert magic here
        return []