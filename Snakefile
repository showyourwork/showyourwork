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


# Return active figure scripts (.py)
def figure_scripts(wildcards):

    # Sleeps until `check_figure_scripts` has run
    checkpoints.check_figure_scripts.get(**wildcards)

    # Get the script: files dict
    with open(Path(".cortex") / "data" / "figure_scripts.json", "r") as f:
        scripts = json.load(f)

    # We need a list of script names
    scripts = scripts.keys()

    # Check the git status of each one
    check_git_status(scripts)

    return scripts


# Return active figure files (.pdf, .png, etc)
def figure_files(wildcards):

    # Sleeps until `check_figure_scripts` has run
    checkpoints.check_figure_scripts.get(**wildcards)

    # Get the script: files dict
    with open(Path(".cortex") / "data" / "figure_scripts.json", "r") as f:
        scripts = json.load(f)

    # We need a list of file names
    files = [file for files in scripts.values() for file in files]

    return files


# Generate figures
rule figures:
    input: "figures/{script}.py"
    output: "figures/{script}.{ext}"
    wildcard_constraints:
        ext="(pdf|png)"
    shell: 
        "cd figures && python {wildcards.script}.py"


# Alternative figure rule (if file names use underscore syntax)
rule figures_alt:
    input: "figures/{script}.py"
    output: "figures/{script}_{suff}.{ext}"
    wildcard_constraints:
        ext="(pdf|png)"
    shell: 
        "cd figures && python {wildcards.script}.py"


# Compile the PDF
rule pdf:
    input: "tex/ms.tex", "tex/bib.bib", figure_files
    output: "ms.pdf"
    run: 
        # Compile the PDF
        if config["debug"]:
            shell("cd tex && tectonic --keep-logs --print ms.tex")
        else:
            shell("cd tex && tectonic --keep-logs ms.tex")
        shell("mv tex/ms.pdf .")


# Remove all temporary files
onsuccess:
    TEMP_FILES = STYLE_FILES_OUT
    TEMP_FILES += glob.glob(str(Path("figures") / "*.py.cortex"))
    for suff in TEMP_SUFFS:
        for f in glob.glob(str(Path("tex") / suff)):
            if "cortex.log" not in f:
                TEMP_FILES.append(f)
    for file in TEMP_FILES:
        os.remove(file)