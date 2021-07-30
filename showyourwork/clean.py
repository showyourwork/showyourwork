from .constants import *
from .utils import glob
import os
import shutil


def clean():
    for ext in FIGURE_EXTENSIONS:
        for file in glob(USER / "figures" / f"*.{ext}"):
            os.remove(file)
    for file in [
        USER / "ms.pdf",
        USER / "dag.pdf",
        USER / ".Snakefile",
        USER / ".helpers.smk",
    ]:
        if file.exists():
            os.remove(file)
    if (USER / ".showyourwork").exists():
        shutil.rmtree(USER / ".showyourwork")


def Clean():
    clean()
    if (USER / ".snakemake").exists():
        shutil.rmtree(USER / ".snakemake")
