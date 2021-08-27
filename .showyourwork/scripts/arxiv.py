import subprocess
import shutil
from pathlib import Path
import os
import tarfile


# Params defined in `../rules/pdf.smk`
verbose = snakemake.params["verbose"]
figexts = snakemake.params["figexts"]
TEMP = snakemake.params["TEMP"]
TEX = snakemake.params["TEX"]
FIGURES = snakemake.params["FIGURES"]
SYWTEXFILE = snakemake.params["SYWTEXFILE"]


# Run tectonic to get the .bbl file
tectonic_args = ["--keep-intermediates"]
if verbose:
    tectonic_args += ["--print"]
else:
    tectonic_args += ["--chatter", "minimal"]
subprocess.check_call(
    ["tectonic"] + tectonic_args + [TEX / "{}.tex".format(SYWTEXFILE)]
)

# Remove all output except the .bbl and .tex files
for file in ["__latexindent_temp.tex", ".showyourwork-ms.*"]:
    for file in TEX.glob(file):
        if file.name not in [".showyourwork-ms.bbl", ".showyourwork-ms.tex"]:
            os.remove(file)

# Copy the `tex` folder over to a temporary location
if (TEMP / "arxiv").exists():
    shutil.rmtree(TEMP / "arxiv")
shutil.copytree(TEX, TEMP / "arxiv")

# Rename our temporary files to `ms.*`
for file in (TEMP / "arxiv").glob(".showyourwork-ms.*"):
    os.rename(file, str(file).replace(".showyourwork-ms", "ms"))

# Remove additional unnecessary files
for file in [
    "ms.pdf",
    "*.bib",
    "sywxml.sty",
    "**/*.py",
    "**/*matplotlibrc",
    "**/*.gitignore",
]:
    for file in (TEMP / "arxiv").glob(file):
        os.remove(file)

# Tar it up
with tarfile.open("arxiv.tar.gz", "w:gz") as tar:
    tar.add(TEMP / "arxiv", arcname=os.path.sep)
shutil.rmtree(TEMP / "arxiv")
