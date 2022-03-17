"""
Miscellaneous functions for interfacing with git.

"""
import subprocess
from pathlib import Path


def get_repo_url():
    """
    Return the repository URL.

    """
    try:
        url = (
            subprocess.check_output(
                ["git", "config", "--get", "remote.origin.url"],
                stderr=subprocess.DEVNULL,
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
    """
    Return the current repository branch name.

    """
    try:
        branch = (
            subprocess.check_output(
                ["git", "branch", "--show-current"], stderr=subprocess.DEVNULL
            )
            .decode()
            .replace("\n", "")
        )
    except Exception as e:
        branch = "unknown"
    return branch


def get_repo_sha():
    """
    Return the SHA for the current git commit.

    """
    try:
        sha = (
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
            )
            .decode()
            .replace("\n", "")
        )
    except Exception as e:
        sha = "unknown"
    return sha


def get_repo_tag():
    """
    Return a tag name if the HEAD corresponds to a tagged version.

    """
    try:
        tag = (
            subprocess.check_output(
                ["git", "describe", "--exact-match", "--tags", "HEAD"],
                stderr=subprocess.DEVNULL
            )
            .decode()
            .strip()
        )
    except Exception:
        tag = ""

    return tag


def get_script_status(script):
    """
    Return an error code corresponding to the git status of a given script.

    """
    # Check if the file exists
    script = Path(script)
    if not script.exists():
        error = ScriptDoesNotExist
    else:
        # Check if the file is version controlled
        try:
            status = (
                subprocess.check_output(
                    ["git", "ls-files", "--error-unmatch", script],
                    cwd=abspaths.user,
                    stderr=subprocess.DEVNULL,
                )
                .decode()
                .replace("\n", "")
            )
        except Exception as e:
            # File is not version controlled
            error = ScriptNotVersionControlled
        else:
            # Check if the file has uncommitted changes
            try:
                status = (
                    subprocess.check_output(
                        ["git", "status", "-s", script],
                        cwd=abspaths.user,
                        stderr=subprocess.DEVNULL,
                    )
                    .decode()
                    .replace("\n", "")
                )
                if len(status):
                    raise Exception("Uncommitted changes!")
            except Exception as e:
                # Uncommited changes
                error = ScriptHasUncommittedChanges
            else:
                # File is good!
                error = ScriptUpToDate
    return error