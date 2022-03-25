from showyourwork import paths
from showyourwork.config import parse_config
from showyourwork.logging import setup_logging, clear_errors
import snakemake
from pathlib import Path
import sys
import os
import jinja2


# Require Snakemake >= this version
snakemake.utils.min_version("6.15.5")


# Working directory is the top level of the user repo
workdir: paths.user().repo.as_posix()


# User config. Allow Jinja2 templating syntax.
with open(paths.user().temp / "showyourwork.yml", "w") as f:
    env = jinja2.Environment(loader=jinja2.FileSystemLoader("."))
    print(env.get_template("showyourwork.yml").render(), file=f)
configfile: (paths.user().temp / "showyourwork.yml").as_posix()


# Report template
report: "report/preprocess.rst"


# Clear errors from past builds
clear_errors()


# Parse the config file
parse_config()


# Set up custom logging
setup_logging(
    verbose=config["verbose"], 
    logfile=paths.user().logs / "preprocess.log"
)


# Hack to make the configfile generation the default rule
rule syw__main:
    input:
        config["config_json"]


# Include all other rules
include: "rules/preprocess.smk"