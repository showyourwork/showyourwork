import subprocess
from pathlib import Path
import shutil
import sys
import os


# The cached output directory
OUTDIR = Path(__file__).resolve().parents[4] / "src/tex/figures"


def get_modified_files(commit="HEAD^"):
    """
    Return all files that changed since `commit`.

    """
    return [
        Path(file)
        for file in (
            subprocess.check_output(
                ["git", "diff", "HEAD", commit, "--name-only"],
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .split("\n")
        )
        if len(file)
    ]


def restore_cache():
    """
    Runs after restoring the cache using @actions/cache.

    """
    # Give everything in the output dir a fresh timestamp
    for file in OUTDIR.rglob("*"):
        file.touch()

    # Get the commit when the files were cached
    try:
        with open(OUTDIR / "last_commit_sha.txt", "r") as f:
            commit = f.readlines()[0].replace("\n", "")
    except FileNotFoundError:
        print("Cache info not found.")

        # DEBUG
        print(OUTDIR)
        print(OUTDIR.rglob("*"))

        return

    # Get all files modified since that commit
    try:
        modified_files = get_modified_files(commit)
    except Exception as e:
        # Can potentially fail if force-pushing is involved!
        print(e)
        return

    # Give all modified files a fresher timestamp than the outputs
    # to trick Snakemake into re-running them.
    for file in modified_files:
        file.touch()


def update_cache():
    """
    Runs before updating the cache using @actions/cache.

    """
    # Store the current commit
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode()
    with open(OUTDIR / "last_commit_sha.txt", "w") as f:
        print(commit, file=f)

    # DEBUG
    print(OUTDIR)
    print(OUTDIR.rglob("*"))


if __name__ == "__main__":
    assert len(sys.argv) == 2, "Incorrect number of args to `cache.py`."
    if sys.argv[1] == "--restore":
        restore_cache()
    elif sys.argv[1] == "--update":
        update_cache()
    else:
        raise ValueError(f"Invalid argument: {sys.argv[1]}")
