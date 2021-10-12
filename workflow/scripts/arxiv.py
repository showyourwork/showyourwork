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
TECTONIC = snakemake.params["TECTONIC"]
arxiv_tarball_exclude = [
    file
    for file in snakemake.params["arxiv_tarball_exclude"].split(",")
    if len(file) > 0
]
arxiv_tarball_exclude += ["sywxml.sty"]


# Run tectonic to get the .bbl file
tectonic_args = ["--keep-intermediates"]
if verbose:
    tectonic_args += ["--print"]
else:
    tectonic_args += ["--chatter", "minimal"]
subprocess.check_call([TECTONIC] + tectonic_args + [TEX / "{}.tex".format(SYWTEXFILE)])

# Remove all output except the .bbl, .tex, and .pdf files
for file in ["__latexindent_temp.tex", ".showyourwork-ms.*"]:
    for file in TEX.glob(file):
        if file.name not in [
            ".showyourwork-ms.bbl",
            ".showyourwork-ms.tex",
            ".showyourwork-ms.pdf",
        ]:
            os.remove(file)

# Copy the `tex` folder over to a temporary location
if (TEMP / "arxiv").exists():
    shutil.rmtree(TEMP / "arxiv")
shutil.copytree(TEX, TEMP / "arxiv")

# Rename our temporary files to `ms.*`
for file in (TEMP / "arxiv").glob(".showyourwork-ms.*"):
    os.rename(file, str(file).replace(".showyourwork-ms", "ms"))

# Remove additional unnecessary files
for file in arxiv_tarball_exclude:
    for file in (TEMP / "arxiv").glob(file):
        try:
            os.remove(file)
        except IsADirectoryError:
            shutil.rmtree(file)

# Tar it up
with tarfile.open("arxiv.tar.gz", "w:gz") as tar:
    tar.add(TEMP / "arxiv", arcname=os.path.sep)
shutil.rmtree(TEMP / "arxiv")
