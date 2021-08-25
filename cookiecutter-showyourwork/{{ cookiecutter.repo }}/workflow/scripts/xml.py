import subprocess


# Params defined in `../rules/tree.smk`
verbose = snakemake.params["verbose"]
TEMP = snakemake.params["TEMP"]
TEX = snakemake.params["TEX"]
TMPTEXFILE = snakemake.params["TMPTEXFILE"]


# Generate the PDF
tectonic_args = ["-r", "0", "-o", TEMP]
if verbose:
    tectonic_args += ["--print"]
else:
    tectonic_args += ["--chatter", "minimal"]
subprocess.check_call(
    ["tectonic"] + tectonic_args + [TEX / "{}.tex".format(TMPTEXFILE)]
)
