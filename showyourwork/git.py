"""
Miscellaneous functions for interfacing with git.

"""
import subprocess


def get_repo_root():
    """
    Return the path to the repository root.

    """
    try:
        root = (
            subprocess.check_output(
                ["git", "rev-parse", "--show-toplevel"],
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .replace("\n", "")
        )
    except Exception as e:
        root = "None"
    return root


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