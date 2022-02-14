import subprocess
from pathlib import Path
import shutil
import sys
import os


# File containing the previous commit SHA
LAST_COMMIT_SHA = Path("src") / "tex" / "figures" / "last_commit_sha.txt"


def get_modified_files(commit="HEAD^"):
    """
    Return all files that changed since `commit`.

    """
    return [
        Path(file).resolve()
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
    # Reset all timestamps across the repo.
    # This ensures no file is newer than any other file, which tricks
    # Snakemake into thinking everything is up to date.
    files = (
        list(Path(".").glob("*"))
        + list((Path(".") / ".showyourwork").rglob("*"))
        + list((Path(".") / "src").rglob("*"))
    )
    for file in files:
        os.utime(file, (0, 0))

    # Get the commit when the files were cached
    try:
        with open(LAST_COMMIT_SHA, "r") as f:
            commit = f.readlines()[0].replace("\n", "")
    except FileNotFoundError:
        print("Cache info not found.")
        return

    # Get all files modified since that commit
    try:
        modified_files = get_modified_files(commit)
    except Exception as e:
        # Can potentially fail if force-pushing is involved!
        print(e)
        return

    # Give all modified files a fresher timestamp than everything else.
    # This will trick Snakemake into re-generating any outputs downstream
    # of these files.
    for file in modified_files:
        print(f"Refreshing timestamp for modified file: {file}")
        file.touch()


def update_cache():
    """
    Runs before updating the cache using @actions/cache.

    TODO: RELATIVE PATHS ARE WRONG. Must be relative to repo root, not `showyourwork`!


    """
    # Store the current commit
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode()
    with open(LAST_COMMIT_SHA, "w") as f:
        print(commit, file=f)


if __name__ == "__main__":
    assert len(sys.argv) == 2, "Incorrect number of args to `cache.py`."
    assert Path(
        "showyourwork/workflow/Snakefile"
    ).exists(), "File `cache.py` must be run from the top level of the repo."
    if sys.argv[1] == "--restore":
        restore_cache()
    elif sys.argv[1] == "--update":
        update_cache()
    else:
        raise ValueError(f"Invalid argument: {sys.argv[1]}")
