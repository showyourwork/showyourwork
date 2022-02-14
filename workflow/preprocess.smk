import snakemake
from pathlib import Path
import sys
import jinja2


# Require Snakemake >= this version
snakemake.utils.min_version("6.7.0")


# Add our utils module to the path
HERE = Path(snakemake.workflow.workflow.current_basedir).absolute()
sys.path.insert(1, str(HERE.parents[0]))
from utils import paths, parse_config, setup_logging, clear_errors


# Working directory is the top level of the user repo
workdir: paths.user.as_posix()



# User config. Allow Jinja2 templating syntax.
with open(paths.temp / "showyourwork.yml", "w") as f:
    env = jinja2.Environment(loader=jinja2.FileSystemLoader("."))
    print(env.get_template("showyourwork.yml").render(), file=f)
configfile: (paths.temp / "showyourwork.yml").as_posix()


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
rule syw__main:
    input:
        config["config_json"]


# Include all other rules
include: "rules/preprocess.smk"