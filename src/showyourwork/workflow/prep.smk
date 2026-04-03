"""
The Snakefile for the article pre-processing step.

"""

from showyourwork import paths, overleaf
from showyourwork.config import render_config, parse_config, get_run_type
from showyourwork.patches import patch_snakemake_logging
from showyourwork.git import get_repo_branch
import os
import shutil


# Working directory is the top level of the user repo
workdir: paths.user().repo.as_posix()


# What kind of run is this? (clean, build, etc.)
run_type = get_run_type()


# User config. Allow Jinja2 templating syntax.
render_config()


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


# Overleaf sync: pull in changes before DAG evaluation so that
# modified tex files are detected by Snakemake.  Previously this lived
# in the ``onstart`` handler, which is skipped when Snakemake decides
# "nothing to be done" – meaning Overleaf edits were silently ignored
# until the user ran ``showyourwork clean`` (see issue #603).
# The env-var sentinel guards against Snakemake evaluating the
# Snakefile more than once per invocation.
if (
    run_type == "preprocess"
    and get_repo_branch() == "main"
    and not os.environ.get("_SYW_OVERLEAF_PULLED")
):
    _overleaf_changed = overleaf.pull_files(
        config["overleaf"]["pull"],
        config["overleaf"]["id"],
        commit_changes=not config["github_actions"],
        push_changes=config["github_actions"]
        and config["overleaf"]["gh_actions_sync"],
    )
    os.environ["_SYW_OVERLEAF_PULLED"] = "1"

    # If Overleaf files changed, remove build outputs so that the
    # subsequent build.smk invocation is forced to recompile the PDF.
    # Without this, Snakemake sees the existing ms.pdf as up-to-date
    # because the Overleaf pull happened outside its tracked DAG.
    if _overleaf_changed:
        _ms_pdf = paths.user().repo / (config.get("ms_name", "ms") + ".pdf")
        if _ms_pdf.exists():
            _ms_pdf.unlink()
        if paths.user().compile.exists():
            shutil.rmtree(paths.user().compile)


# Hack to make the configfile generation the default rule
rule syw__main:
    input:
        config["config_json"],


# Include all other rules
include: "rules/common.smk"
include: "rules/preprocess.smk"
