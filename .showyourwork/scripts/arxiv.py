import subprocess
import shutil
from pathlib import Path
import os
import tarfile


# Params defined in `../rules/pdf.smk`
verbose = snakemake.params["verbose"]
tar = snakemake.params["tar"]
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

# Copy the `tex` folder over to a temporary location
if (TEMP / "arxiv").exists():
    shutil.rmtree(TEMP / "arxiv")
shutil.copytree(TEX, TEMP / "arxiv")

# Rename temporary files to `ms.*`
for file in (TEMP / "arxiv").glob(".showyourwork-ms.*"):
    os.rename(file, str(file).replace(".showyourwork-ms", "ms"))

# Remove garbage / unnecessary files
for file in ["ms.pdf", "__latexindent_temp.tex"]:
    if (TEMP / "arxiv" / file).exists():
        os.remove(TEMP / "arxiv" / file)

# Copy over all figures
for ext in figexts:
    for file in FIGURES.glob("*.{}".format(ext)):
        shutil.copy(file, TEMP / "arxiv")

# Tar it up or keep it as a folder?
if tar:
    with tarfile.open("arxiv.tar.gz", "w:gz") as tar:
        tar.add(TEMP / "arxiv", arcname=os.path.sep)
    shutil.rmtree(TEMP / "arxiv")
else:
    shutil.move(TEMP / "arxiv", "arxiv")
