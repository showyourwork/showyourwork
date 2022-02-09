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
