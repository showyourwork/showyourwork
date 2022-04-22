from showyourwork import paths, overleaf
from showyourwork.config import parse_config, get_run_type
from showyourwork.patches import patch_snakemake_logging
import snakemake
from pathlib import Path
import sys
import os
import jinja2


# Require Snakemake >= this version
snakemake.utils.min_version("6.15.5")


# Working directory is the top level of the user repo
workdir: paths.user().repo.as_posix()


# What kind of run is this? (clean, build, etc.)
run_type = get_run_type()


# User config. Allow Jinja2 templating syntax.
with open(paths.user().temp / "showyourwork.yml", "w") as f:
    env = jinja2.Environment(loader=jinja2.FileSystemLoader("."))
    print(env.get_template("showyourwork.yml").render(), file=f)
configfile: (paths.user().temp / "showyourwork.yml").as_posix()


# Report template
report: "report/preprocess.rst"


# Remove temp flags
for file in paths.user().flags.glob("*"):
    file.unlink()


# Remove old logs
for file in paths.user().logs.glob("*.log"):
    file.unlink()


# Set up custom logging for Snakemake
patch_snakemake_logging()


# Parse the config file
parse_config()


# Hack to make the configfile generation the default rule
rule syw__main:
    input:
        config["config_json"]


# Include all other rules
include: "rules/preprocess.smk"


onstart:

    
    # Overleaf sync: pull in changes
    if run_type == "preprocess":
        overleaf.pull_files(
            config["overleaf"]["pull"], 
            config["overleaf"]["id"], 
            auto_commit=config["overleaf"]["auto-commit"] and config["github_actions"]
        )