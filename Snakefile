import glob
from pathlib import Path
import json
import os
import shutil


# Import cortex
import sys
sys.path.insert(1, ".cortex")
from cortex import gen_pdf, get_user_metadata, get_script_metadata, get_metadata


# Config
configfile: "config.yml"
debug = config["debug"]
figext = "({})".format("|".join(config["figure_extensions"]))

rule user_meta:
    output: ".cortex/data/user.json"
    run: get_user_metadata(clobber=True)

rule scripts_meta:
    input: ".cortex/data/user.json"
    output: ".cortex/data/scripts.json"
    run: get_script_metadata(clobber=True)

rule meta:
    input: ".cortex/data/user.json", ".cortex/data/scripts.json"
    output: ".cortex/data/meta.json"
    run: get_metadata(clobber=True)

rule pdf:
    input: "tex/ms.tex", "tex/bib.bib", ".cortex/data/meta.json"
    output: "ms.pdf"
    run: gen_pdf()