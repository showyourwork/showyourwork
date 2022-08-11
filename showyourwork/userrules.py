import requests

from . import exceptions
from .config import get_run_type
from .git import get_repo_branch
from .logging import get_logger
from .patches import patch_snakemake_cache

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

    # Remind the user to publish the deposit
    if cache_sandbox_doi and get_run_type == "build":
        get_logger().warn(
            f"Zenodo cache not yet published for this repository."
        )

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

    snakemake.workflow.config["cached_deps"] = cached_deps
