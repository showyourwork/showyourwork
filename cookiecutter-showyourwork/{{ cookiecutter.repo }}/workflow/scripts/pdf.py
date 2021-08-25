import subprocess
import shutil


# Params defined in `../rules/pdf.smk`
verbose = snakemake.params["verbose"]
TEMP = snakemake.params["TEMP"]
TEX = snakemake.params["TEX"]
SYWTEXFILE = snakemake.params["SYWTEXFILE"]


# Generate the PDF
tectonic_args = ["-o", TEMP]
if verbose:
    tectonic_args += ["--print"]
else:
    tectonic_args += ["--chatter", "minimal"]
subprocess.check_call(
    ["tectonic"] + tectonic_args + [TEX / "{}.tex".format(SYWTEXFILE)]
)
shutil.copy(TEMP / "{}.pdf".format(SYWTEXFILE), "ms.pdf")
