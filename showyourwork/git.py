"""
Miscellaneous functions for interfacing with git.

"""
from .subproc import get_stdout


def callback(code, stdout, stderr):
    if code != 0:
        return "unknown"
    else:
        return stdout.replace("\n", "")


def get_repo_root():
    """
    Return the path to the repository root.

    """
    return get_stdout(
        ["git", "rev-parse", "--show-toplevel"], callback=callback
    )


def get_commit_message():
    return get_stdout(
        "git log -1 --pretty=%B", shell=True, callback=callback
    )


def get_repo_url():
    """
    Return the repository URL.

    """
    url = get_stdout(
        ["git", "config", "--get", "remote.origin.url"], callback=callback
    )
    if url.endswith(".git"):
        url = url[:-4]
    return url


def get_repo_branch():
    """
    Return the current repository branch name.

    """
    return get_stdout(["git", "branch", "--show-current"], callback=callback)


def get_repo_slug():
    """
    Return the current repository slug ("user/repo").

    """
    return "/".join(get_repo_url().split("/")[-2:])


def get_repo_sha():
    """
    Return the SHA for the current git commit.

    """
    return get_stdout(["git", "rev-parse", "HEAD"], callback=callback)


def get_repo_tag():
    """
    Return a tag name if the HEAD corresponds to a tagged version.

    """
    tag = get_stdout(
        ["git", "describe", "--exact-match", "--tags", "HEAD"],
        callback=callback,
    ).strip()
    if tag == "unknown":
        tag = ""
    return tag