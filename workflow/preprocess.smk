import snakemake
from pathlib import Path
import sys


# Require Snakemake >= this version
snakemake.utils.min_version("6.7.0")


# Add our utils module to the path
HERE = Path(snakemake.workflow.workflow.current_basedir).absolute()
sys.path.insert(1, str(HERE.parents[0]))
from utils import paths, parse_config, setup_logging, clear_errors


# Working directory is the top level of the user repo
workdir: paths.user.as_posix()


# User config
configfile: (paths.user / "showyourwork.yml").as_posix()


# Report template
report: "report/preprocess.rst"


# Clear errors from past builds
clear_errors()


# Parse the config file
parse_config()


# Set up custom logging
setup_logging(
    debug=config["debug"], 
    verbose=config["verbose"], 
    logfile=paths.logs / "preprocess.log"
)


# Hack to make the configfile generation the default rule
rule main:
    input:
        config["config_json"]


# Include all other rules
include: "rules/preprocess.smk"