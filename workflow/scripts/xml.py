"""
Builds the article to get the XML tree.
This script is called from the ``xml`` rule.

"""
import subprocess
import sys
from pathlib import Path


# Hack to import our custom exception
sys.path.insert(0, str(Path(__file__).parents[1] / "rules"))
from exceptions import ShowyourworkException


# Params defined in `../rules/xml.smk`
verbose = snakemake.params["verbose"]
TEMP = snakemake.params["TEMP"]
TEX = snakemake.params["TEX"]
TMPTEXFILE = snakemake.params["TMPTEXFILE"]
TECTONIC = snakemake.params["TECTONIC"]
EXCEPTIONFILE = snakemake.params["EXCEPTIONFILE"]


# Build the LaTeX document to get the XML tree
tectonic_args = ["-r", "0", "-o", TEMP]
if verbose:
    tectonic_args += ["--print"]
else:
    tectonic_args += ["--chatter", "minimal"]
result = subprocess.run(
    [TECTONIC] + tectonic_args + [TEX / "{}.tex".format(TMPTEXFILE)],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
)
if result.returncode > 0:
    raise ShowyourworkException(
        result.stderr.decode("utf-8"),
        exception_file=EXCEPTIONFILE,
        script="xml.py",
        rule_name="xml",
        brief="An error occurred during the initial build of your TeX file.",
        context="This is the initial build step for your article, in which "
        "showyourwork builds an XML tree of your figure environments to "
        "determine the relationship between figure scripts and figure files. "
        "The `src/ms.tex` file is copied to a temporary file called "
        "`.showyourwork-xml-ms.tex`, so you might see references to that file "
        "in the messages below. Note that line numbers may be offset by one "
        "relative to your original `ms.tex` file.",
    )
if verbose:
    print(result.stdout.decode("utf-8"))
    print(result.stderr.decode("utf-8"))

# Add <HTML></HTML> tags to the XML file
try:
    with open(TEMP / "showyourwork.xml", "r") as f:
        contents = f.read()
except FileNotFoundError:
    contents = ""
contents = "<HTML>\n" + contents + "</HTML>"
with open(TEMP / "showyourwork.xml", "w") as f:
    print(contents, file=f)