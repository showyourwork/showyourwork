"""
Creates the output tarball for easy arXiv upload.
This script is called from the ``arxiv`` rule.

"""
import subprocess
import shutil
from pathlib import Path
import os
import tarfile


# Params defined in `../rules/arxiv.smk`
verbose = snakemake.params["verbose"]
figexts = snakemake.params["figexts"]
TEMP = snakemake.params["TEMP"]
SRC = snakemake.params["SRC"]
FIGURES = snakemake.params["FIGURES"]
SYWTEXFILE = snakemake.params["SYWTEXFILE"]
TECTONIC = snakemake.params["TECTONIC"]
ZENODO_FILES = snakemake.params["ZENODO_FILES"]
arxiv_tarball_exclude = snakemake.params["arxiv_tarball_exclude"]


# Paths that are ALWAYS excluded from the tarball
SHOWYOURWORK_EXCLUDE = [
    "src/**/*.py",
    "src/**/*.ipynb",
    "src/**/matplotlibrc",
    "src/**/.gitignore",
    "src/**/__pycache__",
    "src/**/.ipynb_checkpoints",
    "src/**/__latexindent_temp.tex",
    "src/**/sywxml.sty",
    "src/**/*.zenodo",
    "src/**/.DS_Store",
    "src/data",
]


# Run tectonic to get the .bbl file
tectonic_args = ["--keep-intermediates"]
if verbose:
    tectonic_args += ["--print"]
else:
    tectonic_args += ["--chatter", "minimal"]
subprocess.check_call([TECTONIC] + tectonic_args + [SRC / "{}.tex".format(SYWTEXFILE)])


# Remove all output except the .bbl and .tex files
for file in SRC.glob(".showyourwork-ms.*"):
    if file.name not in [".showyourwork-ms.bbl", ".showyourwork-ms.tex"]:
        os.remove(file)


# Files to be excluded
exclude = []
for name in arxiv_tarball_exclude + SHOWYOURWORK_EXCLUDE + ZENODO_FILES:
    for file in Path(".").glob(name):
        if file.is_dir():
            exclude.extend(file.rglob("*"))
        else:
            exclude.append(file)


# Recursively find all files and copy over everything
# (with the above exclusions) to `.showyourwork/arxiv`
if (TEMP / "arxiv").exists():
    shutil.rmtree(TEMP / "arxiv")
for file in Path(SRC).rglob("*"):
    if (not file.is_dir()) and (not file in exclude):
        dest = TEMP / "arxiv" / file.relative_to(SRC)
        os.makedirs(dest.parents[0], exist_ok=True)
        shutil.copy(file, dest)


# Rename our temporary files to `ms.*`
for file in (TEMP / "arxiv").glob(".showyourwork-ms.*"):
    os.rename(file, str(file).replace(".showyourwork-ms", "ms"))


# Tar it up
with tarfile.open("arxiv.tar.gz", "w:gz") as tar:
    tar.add(TEMP / "arxiv", arcname="arxiv")
shutil.rmtree(TEMP / "arxiv")
