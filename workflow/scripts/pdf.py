"""
Builds the article PDF.
This script is called from the ``pdf`` rule.

"""
import subprocess
import shutil
import sys
from pathlib import Path


# Hack to import our custom exception
sys.path.insert(0, str(Path(__file__).parents[1] / "helpers"))
from helpers.exceptions import ShowyourworkException


# Params defined in `../rules/pdf.smk`
verbose = snakemake.params["verbose"]
TEMP = snakemake.params["TEMP"]
TEX = snakemake.params["TEX"]
SYWTEXFILE = snakemake.params["SYWTEXFILE"]
TECTONIC = snakemake.params["TECTONIC"]
EXCEPTIONFILE = snakemake.params["EXCEPTIONFILE"]

# Generate the PDF
tectonic_args = ["-o", TEMP]
if verbose:
    tectonic_args += ["--print"]
else:
    tectonic_args += ["--chatter", "minimal"]
result = subprocess.run(
    [TECTONIC] + tectonic_args + [TEX / "{}.tex".format(SYWTEXFILE)],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
)
if result.returncode > 0:
    raise ShowyourworkException(
        result.stderr.decode("utf-8"),
        exception_file=EXCEPTIONFILE,
        script="pdf.py",
        rule_name="pdf",
        brief="An error occurred during the final build of your TeX file.",
        context="This is the final build step for your article that "
        "generates your article PDF `ms.pdf`. "
        "The `src/ms.tex` file is copied to a temporary file called "
        "`.showyourwork-ms.tex`, so you might see references to that file "
        "in the messages below. Note that line numbers may be offset by one "
        "relative to your original `ms.tex` file.",
    )
if verbose:
    print(result.stdout.decode("utf-8"))
    print(result.stderr.decode("utf-8"))
shutil.move(TEMP / "{}.pdf".format(SYWTEXFILE), "ms.pdf")