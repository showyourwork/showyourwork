"""
Builds the article to get the XML tree.
This script is called from the ``xml`` rule.

"""
import subprocess


# Params defined in `../rules/xml.smk`
verbose = snakemake.params["verbose"]
TEMP = snakemake.params["TEMP"]
TEX = snakemake.params["TEX"]
TMPTEXFILE = snakemake.params["TMPTEXFILE"]
TECTONIC = snakemake.params["TECTONIC"]


# Build the LaTeX document to get the XML tree
tectonic_args = ["-r", "0", "-o", TEMP]
if verbose:
    tectonic_args += ["--print"]
else:
    tectonic_args += ["--chatter", "minimal"]
subprocess.check_call([TECTONIC] + tectonic_args + [TEX / "{}.tex".format(TMPTEXFILE)])

# Add <HTML></HTML> tags to the XML file
try:
    with open(TEMP / "showyourwork.xml", "r") as f:
        contents = f.read()
except FileNotFoundError:
    contents = ""
contents = "<HTML>\n" + contents + "</HTML>"
with open(TEMP / "showyourwork.xml", "w") as f:
    print(contents, file=f)