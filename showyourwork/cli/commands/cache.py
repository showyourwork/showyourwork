from ... import paths
from ...subproc import get_stdout
from pathlib import Path
import os
import time


def get_modified_files(commit="HEAD^"):
    """
    Return all files that changed since `commit`.

    """
    return [
        Path(file).resolve()
        for file in (
            get_stdout(["git", "diff", "HEAD", commit, "--name-only"]).split(
                "\n"
            )
        )
        if len(file)
    ]


def cache_restore():
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
    timestamp = time.time()
    for file in files:
        os.utime(file, (timestamp, timestamp))

    # Get the commit when the files were cached
    try:
        with open(paths.user().figures / "last_commit_sha.txt", "r") as f:
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

    # Give all modified files a newer timestamp than everything else.
    # This will trick Snakemake into re-generating any outputs downstream
    # of these files.
    time.sleep(1)
    for file in modified_files:
        print(f"Refreshing timestamp for modified file: {file}")
        file.touch()


def cache_update():
    """
    Runs before updating the cache using @actions/cache.

    """
    # Store the current commit
    commit = get_stdout(["git", "rev-parse", "HEAD"])
    with open(paths.user().figures / "last_commit_sha.txt", "w") as f:
        print(commit, file=f)
