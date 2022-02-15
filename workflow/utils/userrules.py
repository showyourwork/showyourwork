from distutils.command.upload import upload
from . import exceptions
from .logging import get_logger
from .zenodo import download_file_from_draft, upload_file_to_draft
from .config import get_snakemake_variable
import types
import snakemake

__all__ = ["process_user_rules"]


def patch_snakemake_cache():
    """
    Patches the Snakemake cache functionality to

        - Add custom logging messages
        - Attempt to download the cache file from Zenodo on `fetch()`
        - Uploads the cache file to Zenodo on `store()`
        - Prevent the contents of the Snakefile defining a rule from getting
          ingested into the hash for that rule. Showyourwork automatically adds
          the Snakefiles as inputs to user-generated rules to force re-running
          them when their definition changes. However, if a Snakefile defines
          multiple rules (and the hash for the rule depends on the contents of
          the Snakefile), editing the parameters for one rule will force
          re-running _all_ rules, since they'll all have different hashes. To
          circumvent this, we intercept the `hash_file()` function and filter
          out the Snakefile associated with the current rule. Note that changing
          the parameters for the rule _still_ will force re-evaluation, since
          the parameters themselves are also ingested when computing the rule
          hash.

    """
    # Get the showyourwork logger
    logger = get_logger()

    # The instance we'll patch
    output_file_cache = snakemake.workflow.workflow.output_file_cache

    # Make a copy of the original methods
    _fetch = output_file_cache.fetch
    _store = output_file_cache.store
    _hash_file = snakemake.caching.hash.hash_file

    # Define the patches
    def fetch(self, job):
        # If the cache file is a directory, we must tar it up
        # Recall that cacheable jobs can only have a _single_ output
        # (unless using multiext), so checking the first output should suffice
        if job.output[0].is_directory:
            tarball = True
        else:
            tarball = False
        for outputfile, cachefile in self.get_outputfiles_and_cachefiles(job):
            if not cachefile.exists():
                # Attempt to download from Zenodo
                try:
                    logger.info(f"Searching Zenodo file cache: {outputfile}...")
                    download_file_from_draft(cachefile, job.rule.name, tarball=tarball)
                    logger.info(f"Restoring from Zenodo cache: {outputfile}...")
                except exceptions.FileNotFoundOnZenodo:
                    # Cache miss; not fatal
                    logger.warn(f"File not found in cache: {outputfile}.")
            else:
                logger.info(f"Restoring from local cache: {outputfile}...")
                # Always ensure the cached file on Zenodo is the latest
                # local cache hit. (If it is, this is a no-op)
                logger.info(f"Syncing file with Zenodo cache: {outputfile}...")
                upload_file_to_draft(cachefile, job.rule.name, tarball=tarball)

        # Call the original method
        return _fetch(job)

    def store(self, job):
        # Call the original method
        result = _store(job)

        # See note in `fetch()` about tarballs
        if job.output[0].is_directory:
            tarball = True
        else:
            tarball = False

        for outputfile, cachefile in self.get_outputfiles_and_cachefiles(job):
            logger.info(f"Caching output file on Zenodo: {outputfile}...")
            upload_file_to_draft(cachefile, job.rule.name, tarball=tarball)

        return result

    def hash_file(f):
        # Don't hash the Snakefile in which the rule is declared
        job = get_snakemake_variable("job")
        if f == job.rule.snakefile:
            return ""
        else:
            return _hash_file(f)

    # Apply them
    output_file_cache.fetch = types.MethodType(fetch, output_file_cache)
    output_file_cache.store = types.MethodType(store, output_file_cache)
    snakemake.caching.hash.hash_file = hash_file


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
    if snakemake.workflow.config["zenodo_cache"]:
        patch_snakemake_cache()

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
        # and add the env file as an explicit input
        if ur.conda_env:
            if hasattr(ur.conda_env, "is_file") and ur.conda_env.is_file:
                ur.set_input(str(ur.conda_env.file))
            else:
                ur.set_input(ur.conda_env)
        else:
            # All users rules should run in conda envs!
            # TODO
            raise exceptions.MissingCondaEnvironmentInUserRule()

        # Add the Snakefile or .smk file defining the rule as an explicit
        # input (user can modularize their Snakefile for more fine-grained
        # control here)
        ur.set_input(ur.snakefile)

        """
        # Make all user rules cacheable
        # DEPRECATED: User must now explicitly set `cache: True` 
        # in rules whose output is to be cached
        ur.ruleinfo.cache = True
        snakemake.workflow.workflow.cache_rules.add(ur.name)
        try:
            ur.check_caching()
        except snakemake.exceptions.RuleException:
            ur.ruleinfo.cache = False
            snakemake.workflow.workflow.cache_rules.remove(ur.name)
        """