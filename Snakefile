import glob
from pathlib import Path
import json
import os
import shutil
import sys
sys.path.insert(1, ".cortex")
import cortex


# Config
configfile: "config.yml"
debug = config["debug"]
figext = "({})".format("|".join(config["figure_extensions"]))


rule user_meta:
    output: ".cortex/data/user.json"
    run: cortex.get_user_metadata(clobber=True)


checkpoint check_scripts:
    input: "tex/ms.tex"
    output: ".cortex/data/scripts.json"
    run: cortex.get_script_metadata(clobber=True)


def scripts(wildcards):
    # Sleeps until `check_scripts` has run
    checkpoints.check_scripts.get(**wildcards)
    return cortex.get_scripts()


rule meta:
    input: ".cortex/data/user.json", scripts
    output: ".cortex/data/meta.json"
    run: cortex.get_metadata(clobber=True)


rule pdf:
    input: "tex/ms.tex", "tex/bib.bib", ".cortex/data/meta.json"
    output: "ms.pdf"
    run: cortex.gen_pdf()