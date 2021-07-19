from ..constants import *
from ..utils import save_json
import subprocess


def get_repo_metadata(clobber=True):
    """
    Return repository metadata.

    """
    if clobber or not (TEMP / "user" / "repo.json").exists():

        repo = {}
        save = True

        # Repo url
        try:
            repo["url"] = (
                subprocess.check_output(
                    ["git", "config", "--get", "remote.origin.url"],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    cwd=USER,
                )
                .decode()
                .replace("\n", "")
            )
            if repo["url"].endswith(".git"):
                repo["url"] = repo["url"][:-4]
        except Exception as e:
            repo["url"] = "unknown"
            save = False

        # Current branch
        try:
            repo["branch"] = (
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
            repo["branch"] = "unknown"
            save = False

        # Save
        if save:
            save_json(repo, TEMP / PROJECT / "repo.json")

    else:

        # Load from profile
        with open(TEMP / PROJECT / "repo.json", "r") as f:
            repo = json.load(f)

    # Get the current git hash (updated every time)
    try:
        repo["sha"] = (
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
        repo["sha"] = "unknown"

    return repo
