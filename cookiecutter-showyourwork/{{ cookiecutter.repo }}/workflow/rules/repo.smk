import subprocess
import json


def get_repo_url():
    try:
        url = (
            subprocess.check_output(
                ["git", "config", "--get", "remote.origin.url"],
                stderr=subprocess.DEVNULL
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
                stderr=subprocess.DEVNULL
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
                stderr=subprocess.DEVNULL
            )
            .decode()
            .replace("\n", "")
        )
    except Exception as e:
        sha = "unknown"
    return sha


rule repo:
    message:
        "Generating repo metadata..."
    output:
        temp(POSIX(TEMP / "repo.json"))
    priority:
        99
    run:
        repo = {}
        repo["url"] = get_repo_url()
        repo["branch"] = get_repo_branch()
        with open(TEMP / "repo.json", "w") as f:
            print(json.dumps(repo, indent=4), file=f)
