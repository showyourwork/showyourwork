import subprocess
import shutil


# Params defined in `../rules/pdf.smk`
verbose = snakemake.params["verbose"]
TEMP = snakemake.params["TEMP"]
TEX = snakemake.params["TEX"]
SYWTEXFILE = snakemake.params["SYWTEXFILE"]
TECTONIC = snakemake.params["TECTONIC"]

# Generate the PDF
tectonic_args = ["-o", TEMP]
if verbose:
    tectonic_args += ["--print"]
else:
    tectonic_args += ["--chatter", "minimal"]
subprocess.check_call([TECTONIC] + tectonic_args + [TEX / "{}.tex".format(SYWTEXFILE)])
shutil.move(TEMP / "{}.pdf".format(SYWTEXFILE), "ms.pdf")