from .constants import *
from pathlib import Path
import shutil
import subprocess
import os


def dag():
    files = list((ROOT / "workflow").glob("*"))  # NOTE: this includes dotfiles
    for file in files:
        shutil.copy(file, USER)
    dag = subprocess.check_output(
        ["snakemake", "--dag", "--snakefile", ".Snakefile"], cwd=USER
    ).decode()
    with open(USER / "dag.dag", "w") as f:
        print(dag, file=f)
    with open(USER / "dag.pdf", "wb") as f:
        f.write(
            subprocess.check_output(["dot", "-Tpdf", str(USER / "dag.dag")])
        )
    os.remove(USER / "dag.dag")
    for file in files:
        os.remove(USER / Path(file).name)
