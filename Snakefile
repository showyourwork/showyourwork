import glob
from pathlib import Path
import json
import os
import shutil
import sys
sys.path.insert(1, ".cortex")
from generate_sty import generate_sty
from check_git_status import check_git_status
from infer_figure_scripts import infer_figure_scripts
from user_prompt import user_prompt


# Config
configfile: "config.yml"
debug = config["debug"]
tectonic_flags = []
if debug:
    tectonic_flags += ["--print"]
if not config["quiet"]:
    tectonic_flags += ["--keep-logs"] 
tectonic_flags = " ".join(tectonic_flags)
figext = "({})".format("|".join(config["figure_extensions"]))


# Generate the metadata file from direct user input
if not (Path(".cortex") / "data" / "meta.json").exists():
    user_prompt()


# Generate `cortex.sty` from the template
generate_sty()


# Temporary LaTeX files
TEMP_SUFFS = ["*.log", "*.blg"]
STYLE_FILES_IN = glob.glob(str(Path(".cortex") / "styles" / "*"))
STYLE_FILES_OUT = [
    str(Path("tex") / Path(file).name) for file in STYLE_FILES_IN
]
for f in STYLE_FILES_IN:
    shutil.copy(f, "tex/")


# Infer the active figure scripts whenever `ms.tex` changes
checkpoint check_figure_scripts:
    input: 
        "tex/ms.tex"
    output: 
        ".cortex/data/figure_scripts.json"
    run: 
        infer_figure_scripts()


# Return active figure script metadata (.py.cortex)
def figure_metadata(wildcards):

    # Sleeps until `check_figure_scripts` has run
    checkpoints.check_figure_scripts.get(**wildcards)

    # Get the script: files dict
    with open(Path(".cortex") / "data" / "figure_scripts.json", "r") as f:
        scripts = json.load(f)

    # We need a list of script names
    scripts = scripts.keys()

    # Check the git status of each one
    check_git_status(scripts)

    # Return the metadata files
    return [
        str(Path(".cortex") / "data" / "{}.cortex".format(Path(script).name)) 
        for script in scripts
    ]


# Return active figure files (.pdf, .png, etc)
def figure_files(wildcards):

    # Sleeps until `check_figure_scripts` has run
    checkpoints.check_figure_scripts.get(**wildcards)

    # Get the script: files dict
    with open(Path(".cortex") / "data" / "figure_scripts.json", "r") as f:
        scripts = json.load(f)

    # Return all figure files
    return [file for files in scripts.values() for file in files]


# Generate figures
rule figures:
    input: "figures/{script}.py"
    output: "figures/{script}.{figext}"
    wildcard_constraints:
        figext=figext
    shell: 
        "cd figures && python {wildcards.script}.py"


# Alternative figure rule (if file names use underscore syntax)
rule figures_alt:
    input: "figures/{script}.py"
    output: "figures/{script}_{suff}.{figext}"
    wildcard_constraints:
        figext=figext
    shell: 
        "cd figures && python {wildcards.script}.py"


# Compile the PDF
rule pdf:
    input: "tex/ms.tex", "tex/bib.bib", figure_files, figure_metadata
    output: "ms.pdf"
    shell: f"cd tex && tectonic {tectonic_flags} ms.tex && mv ms.pdf ../"


# Remove all temporary files
onsuccess:
    if not debug:
        TEMP_FILES = STYLE_FILES_OUT
        TEMP_FILES += glob.glob(str(Path("figures") / "*.py.cortex"))
        for suff in TEMP_SUFFS:
            for f in glob.glob(str(Path("tex") / suff)):
                if "cortex.log" not in f:
                    TEMP_FILES.append(f)
        for file in TEMP_FILES:
            os.remove(file)