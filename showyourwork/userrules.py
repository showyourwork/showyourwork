from . import exceptions
from .logging import get_logger
from .patches import patch_snakemake_cache
from .config import get_run_type
from .git import get_repo_branch
import json
import requests

try:
    import snakemake
except ModuleNotFoundError:
    snakemake = None


__all__ = ["process_user_rules"]


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
    cached_deps = []
    branch = get_repo_branch()
    if snakemake.workflow.config["cache"][branch]:
        patch_snakemake_cache(snakemake.workflow.config["cache"][branch])

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
            ur.set_input(ur.script)
        elif ur.notebook:
            ur.set_input(ur.notebook)
        elif ur.is_run:
            raise exceptions.RunDirectiveNotAllowedInUserRules(ur.name)

        # Record any cached output
        if snakemake.workflow.config["cache"][branch]:
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
    try:

        # Grab the Zenodo id for this repo
        concept_id = snakemake.workflow.config["cache"][branch]

        # Check that the record has been published
        if concept_id:
            zenodo_cache_url = f"https://zenodo.org/record/{concept_id}"
            r = requests.get(f"https://zenodo.org/api/record/{concept_id}")
            data = r.json()
            if r.status_code > 204:
                zenodo_cache_url = None
                get_logger().warn(
                    f"Zenodo cache record {concept_id} "
                    "not yet published for this repository."
                )
                raise Exception()
        else:
            raise Exception()

    except:

        # Fail silently
        pass

    else:

        # Add the metadata to the config (accessed in `pdf.py`)
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