from .constants import *
from .utils import glob
from .clean import clean, Clean
from .cache import restore_cache, update_cache
from .new import new
from .dag import dag
import subprocess
import shutil
import os
from pathlib import Path
import argparse
import sys
import re


def main():
    # Parse command line args
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--clean", action="store_true")
    parser.add_argument("-C", "--Clean", action="store_true")
    parser.add_argument("-n", "--new", action="store_true")
    parser.add_argument("-d", "--dag", action="store_true")

    # Internal arguments
    parser.add_argument(
        "--restore-cache", action="store_true", help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--update-cache", action="store_true", help=argparse.SUPPRESS
    )
    args, snakemake_args = parser.parse_known_args(sys.argv[1:])

    # Subprograms
    cmds = ["clean", "Clean", "new", "dag", "restore_cache", "update_cache"]
    assert (
        sum([int(getattr(args, cmd, False)) for cmd in cmds]) < 2
    ), "Options conflict!"  # Only one is allowed at a time!
    for cmd in cmds:
        if getattr(args, cmd, False):
            exec(f"{cmd}()")
            return

    # Process Snakemake defaults
    cores_set = False
    conda_set = False
    for i, arg in enumerate(snakemake_args):
        if re.match("^-c([0-9]*?)$", arg):
            cores_set = True
        if arg == "--use-conda":
            conda_set = True
        if (arg == "-s") or (arg == "--snakefile"):
            raise ValueError(
                "Arguments `-s` or `--snakefile` are not allowed."
            )

        if arg == "--verbose":
            # TODO: Implement me!
            pass

    if not cores_set:
        snakemake_args.append("-c1")
    if not conda_set:
        snakemake_args.append("--use-conda")
    snakemake_args.extend(["--snakefile", ".Snakefile"])

    # Copy workflow to user directory
    files = list((ROOT / "workflow").glob("*"))  # NOTE: this includes dotfiles
    for file in files:
        shutil.copy(file, USER)

    # Run Snakemake
    subprocess.check_call(["snakemake"] + snakemake_args, cwd=USER)

    # Remove tempfiles
    for file in files:
        os.remove(USER / Path(file).name)
    for file in glob(USER / "tex" / "*latexindent*.tex"):
        os.remove(file)
