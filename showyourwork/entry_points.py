from .constants import *
from .utils import glob
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
    parser.add_argument("--clean", action="store_true")
    args, snakemake_args = parser.parse_known_args(sys.argv[1:])

    # Clean
    if args.clean:
        for ext in FIGURE_EXTENSIONS:
            for file in glob(USER / "figures" / f"*.{ext}"):
                os.remove(file)
        if (USER / "ms.pdf").exists():
            os.remove(USER / "ms.pdf")
        return

    # Process Snakemake defaults
    cores_set = False
    conda_set = False
    for arg in snakemake_args:
        if re.match("^-c([0-9]*?)$", arg):
            cores_set = True
        if arg == "--use-conda":
            conda_set = True
    if not cores_set:
        snakemake_args.append("-c1")
    if not conda_set:
        snakemake_args.append("--use-conda")

    # Copy workflow to user directory
    files = glob(ROOT / "workflow" / "*")
    for file in files:
        shutil.copy(file, USER)

    # Run Snakemake
    subprocess.check_output(["snakemake"] + snakemake_args, cwd=USER)

    # Remove tempfiles
    for file in files:
        os.remove(USER / Path(file).name)
    for file in glob(USER / "tex" / "*latexindent*.tex"):
        os.remove(file)
