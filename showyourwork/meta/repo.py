from ..constants import *
from ..utils import save_json
import subprocess
import json


__all__ = ["get_repo_metadata"]


def get_repo_url():
    try:
        url = (
            subprocess.check_output(
                ["git", "config", "--get", "remote.origin.url"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                cwd=USER,
            )
            .decode()
            .replace("\n", "")
        )
        if url.endswith(".git"):
            url = url[:-4]
    except Exception as e:
        url = "unknown"
    return url


def get_repo_branch():
    try:
        branch = (
            subprocess.check_output(
                ["git", "branch", "--show-current"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                cwd=USER,
            )
            .decode()
            .replace("\n", "")
        )
    except Exception as e:
        branch = "unknown"
    return branch


def get_repo_sha():
    try:
        sha = (
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                cwd=USER,
            )
            .decode()
            .replace("\n", "")
        )
    except Exception as e:
        sha = "unknown"
    return sha


def get_repo_metadata(clobber=True):
    """
    Return repository metadata.

    """
    if clobber or not (TEMP / "user" / "repo.json").exists():

        # Get git info
        repo = {}
        repo["url"] = get_repo_url()
        repo["branch"] = get_repo_branch()

        # Save
        save_json(repo, TEMP / PROJECT / "repo.json")

    else:

        # Load from profile
        with open(TEMP / PROJECT / "repo.json", "r") as f:
            repo = json.load(f)

        # Check if the cached info is good
        if repo["url"] == "unknown":
            repo["url"] = get_repo_url()
        if repo["branch"] == "unknown":
            repo["branch"] = get_repo_branch()

    # Get the current git hash (updated every time)
    repo["sha"] = get_repo_sha()

    return repo
