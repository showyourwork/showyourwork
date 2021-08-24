from .constants import *
from .utils import glob, check_repo
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
    parser = argparse.ArgumentParser(
        description="Run this command in the top level of a scientific "
        "article repository "
        "to execute the Snakemake workflow and build the article PDF. "
        "See https://github.com/rodluger/showyourwork for more details. "
        "In addition to the options below, users may provide any of the "
        "arguments accepted by `Snakemake`, as well as a `-v, --verbose` "
        "option to increase the verbosity of both `Snakemake` and `tectonic`. "
    )
    parser.add_argument(
        "-c",
        "--clean",
        action="store_true",
        help="remove all temporary files and workflow outputs",
    )
    parser.add_argument(
        "-C",
        "--Clean",
        action="store_true",
        help="runs `clean` and removes the Snakemake cache",
    )
    parser.add_argument(
        "-n",
        "--new",
        nargs=1,
        metavar="slug",
        help="create a new article repository",
    )
    parser.add_argument(
        "-d",
        "--dag",
        action="store_true",
        help="generate a directed acyclic graph (DAG) of the workflow, "
        "saved as `dag.pdf`",
    )

    # Internal arguments
    parser.add_argument(
        "--restore-cache", action="store_true", help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--update-cache", action="store_true", help=argparse.SUPPRESS
    )
    args, snakemake_args = parser.parse_known_args(sys.argv[1:])

    # Subprograms
    cmds = ["clean", "Clean", "dag", "restore_cache", "update_cache"]
    assert (
        sum([int(not getattr(args, cmd, False) is False) for cmd in cmds]) < 2
    ), "Options conflict!"  # Only one is allowed at a time!
    for cmd in cmds:
        if getattr(args, cmd, False):
            exec(f"{cmd}()")
            return
    if args.new:
        new(args.new[0])
        return

    # Check that the `cwd` is a valid article repo
    check_repo()

    # Process Snakemake defaults
    cores_set = False
    conda_set = False
    for i, arg in enumerate(snakemake_args):
        if re.match("^-c([0-9]*?)$", arg):
            cores_set = True
        if arg == "--use-conda":
            conda_set = True
        if arg == "--verbose":
            snakemake_args.extend(["--config", "verbose=true"])

    if not cores_set:
        snakemake_args.append("-c1")
    if not conda_set:
        snakemake_args.append("--use-conda")

    # Run Snakemake
    subprocess.check_call(["snakemake"] + snakemake_args, cwd=USER)

    # Remove tempfiles
    for file in glob(USER / "tex" / "*latexindent*.tex"):
        os.remove(file)
