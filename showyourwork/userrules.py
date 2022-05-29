from . import exceptions
from .logging import get_logger
from .patches import patch_snakemake_cache, get_snakemake_variable
from .config import get_run_type
from .git import get_repo_branch
import json
import requests

try:
    import snakemake
except ModuleNotFoundError:
    snakemake = None


def process_user_rules():
    """
    Process user-defined Snakemake rules.

    """
    # Get all showyourwork and user rules
    syw_rules = []
    user_rules = []
    for r in snakemake.workflow.workflow.rules:
        if r.name.startswith("syw__"):
            syw_rules.append(r)
        else:
            user_rules.append(r)

    # Patch the Snakemake caching functionality so we
    # can cache things on Zenodo
    branch = get_repo_branch()
    cache_zenodo_doi = snakemake.workflow.config["cache"][branch]["zenodo"]
    cache_sandbox_doi = snakemake.workflow.config["cache"][branch]["sandbox"]
    cached_deps = []
    if cache_zenodo_doi or cache_sandbox_doi:
        patch_snakemake_cache(cache_zenodo_doi, cache_sandbox_doi)

    # Process each user rule
    for ur in user_rules:

        # Set its order > all showyourwork rules
        for sr in syw_rules:
            snakemake.workflow.workflow.ruleorder(ur.name, sr.name)

        # Add a message if it doesn't have one
        if not ur.message:
            ur.message = f"Running user rule {ur.name}..."

        # Add script as an explicit input
        if ur.script:
            script = ur.script
            script = script.replace("{wildcards.", "{")
            ur.set_input(script)
        elif ur.notebook:
            notebook = ur.notebook
            notebook = notebook.replace("{wildcards.", "{")
            ur.set_input(notebook)
        elif ur.is_run:
            raise exceptions.RunDirectiveNotAllowedInUserRules(ur.name)

        # Record any cached output
        if cache_zenodo_doi or cache_sandbox_doi:
            if ur.ruleinfo.cache:
                for file in ur.output:
                    cached_deps.append(str(file))

        # Ensure we're running in a conda env
        if not ur.conda_env:
            ur.conda_env = "environment.yml"

    # No need to go any further if this isn't the main build
    if get_run_type() != "main":
        return

    # Add some metadata containing the link to the Zenodo cache record
    # which we can access on the LaTeX side to provide margin links
    # for figures that depend on cached datasets
    if cache_zenodo_doi:

        # Add the metadata to the config (accessed in `pdf.py`)
        zenodo_cache_url = f"https://zenodo.org/record/{cache_zenodo_doi}"
        labels = {}
        for figscript in snakemake.workflow.config["dependencies"]:
            for dep in snakemake.workflow.config["dependencies"][figscript]:
                if dep in cached_deps:
                    for label in snakemake.workflow.config["labels"]:
                        if (
                            label.endswith("_script")
                            and snakemake.workflow.config["labels"][label]
                            == figscript
                        ):
                            labels[label[:-6] + "cache"] = zenodo_cache_url
        snakemake.workflow.config["labels"].update(labels)

    else:

        if cache_sandbox_doi:
            # There's a Sandbox deposit but no Zenodo deposit; remind
            # the user to publish the deposit
            get_logger().warn(
                f"Zenodo cache not yet published for this repository."
            )
