"""
Miscellaneous functions for interfacing with git.

"""
from .subproc import get_stdout


def callback(code, stdout, stderr):
    if code != 0:
        return "None"
    else:
        return stdout.replace("\n", "")


def get_repo_root():
    """
    Return the path to the repository root.

    """
    return get_stdout(
        ["git", "rev-parse", "--show-toplevel"], callback=callback
    )


def get_repo_url():
    """
    Return the repository URL.

    """
    try:
        url = get_stdout(
            ["git", "config", "--get", "remote.origin.url"], callback=callback
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
        branch = get_stdout(
            ["git", "branch", "--show-current"], callback=callback
        )
    except Exception as e:
        branch = "unknown"
    return branch


def get_repo_sha():
    """
    Return the SHA for the current git commit.

    """
    try:
        sha = get_stdout(["git", "rev-parse", "HEAD"], callback=callback)
    except Exception as e:
        sha = "unknown"
    return sha