from click import exceptions
from . import zenodo
import snakemake

__all__ = ["process_user_rules"]


def process_user_rules():
    """
    Process user-defined Snakemake rules.

    """

    # Get figure, zenodo, and user rules
    syw_rules = []
    user_rules = []
    for r in snakemake.workflow.workflow.rules:
        if r.name.startswith("syw__"):
            syw_rules.append(r)
        else:
            user_rules.append(r)

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
            # TODO
            raise exceptions.RunDirectiveNotAllowedInUserRules()

        # Ensure we're running in a conda env
        if ur.conda_env:
            ur.set_input(ur.conda_env)
        else:
            # All users rules should run in conda envs!
            # TODO
            raise exceptions.MissingCondaEnvironmentInUserRule()

        # Add the user Snakefile as an input
        ur.set_input(ur.snakefile)

        # Make the rule cacheable
        ur.ruleinfo.cache = True
        snakemake.workflow.workflow.cache_rules.add(ur.name)
        try:
            ur.check_caching()
        except snakemake.exceptions.RuleException:
            ur.ruleinfo.cache = False
            snakemake.workflow.workflow.cache_rules.remove(ur.name)
