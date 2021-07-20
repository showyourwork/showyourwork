from .constants import *
import subprocess
from pathlib import Path
from glob import glob
import os
import shutil
import argparse


# Jan 1 2000 00:00:00
ANCIENT = "200001010000.00"


def get_modified_files(commit="HEAD^"):
    """
    Return all files that changed since `commit`.

    """
    return [
        file
        for file in (
            subprocess.check_output(
                ["git", "diff", "HEAD", commit, "--name-only"]
            )
            .decode()
            .split("\n")
        )
        if len(file)
    ]


def restore_cache():
    """
    Restore cached output files.

    """
    # Copy the cached figure files
    for ext in FIGURE_EXTENSIONS:
        for file in glob(
            str(USER / ".cache" / "figures" / "*.{}".format(ext))
        ):
            shutil.copy2(file, USER / "figures")
            print(
                "Restored file figures/{} from cache.".format(Path(file).name)
            )

    # Get the commit when the files were cached
    try:
        with open(USER / ".cache" / "commit", "r") as f:
            commit = f.readlines()[0].replace("\n", "")
    except FileNotFoundError:
        return

    # Get all files modified since that commit
    files = get_modified_files(commit)

    # If an input file has not changed since the commit at which its
    # output was cached, mark it as ancient to prevent Snakemake
    # from re-running it.
    input_scripts = glob(str(USER / "figures" / "*.py"))
    input_scripts += glob(str(USER / "tex" / "*.bib"))
    input_scripts += [str(USER / "tex" / "ms.tex")]
    input_scripts += [str(USER / "environment.yml")]
    for script in input_scripts:
        script_rel = str(Path("figures") / Path(script).name)
        if script_rel not in files:
            subprocess.check_output(
                ["touch", "-a", "-m", "-t", ANCIENT, script]
            )
            print("File `{}` marked as ANCIENT.".format(script_rel))


def update_cache():
    """
    Cache output files.

    """
    # Clear the existing cache
    if (USER / ".cache").exists():
        shutil.rmtree(USER / ".cache")
    os.mkdir(USER / ".cache")
    os.mkdir(USER / ".cache" / "figures")

    # Cache figure files
    for ext in FIGURE_EXTENSIONS:
        for file in glob(str(USER / "figures" / "*.{}".format(ext))):
            shutil.copy2(file, USER / ".cache" / "figures")
            print("Cached file figures/{}.".format(Path(file).name))

    # Store the current commit
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode()
    with open(USER / ".cache" / "commit", "w") as f:
        print(commit, file=f)
