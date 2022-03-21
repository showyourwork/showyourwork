from distutils.command.upload import upload
from . import exceptions
from .logging import get_logger
from .zenodo import download_file_from_zenodo, upload_file_to_zenodo
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

    """
    # Get the showyourwork logger
    logger = get_logger()

    # The instance we'll patch
    output_file_cache = snakemake.workflow.workflow.output_file_cache

    # Make a copy of the original methods
    _fetch = output_file_cache.fetch
    _store = output_file_cache.store

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
                    download_file_from_zenodo(cachefile, job.rule.name, tarball=tarball)
                    logger.info(f"Restoring from Zenodo cache: {outputfile}...")
                except exceptions.FileNotFoundOnZenodo:
                    # Cache miss; not fatal
                    logger.warn(f"File not found in cache: {outputfile}.")
                except exceptions.ZenodoException as e:
                    # NOTE: we treate all Zenodo caching errors as non-fatal
                    logger.error(str(e))
            else:
                logger.info(f"Restoring from local cache: {outputfile}...")
                # Always ensure the cached file on Zenodo is the latest
                # local cache hit. (If it is, this is a no-op)
                logger.info(f"Syncing file with Zenodo cache: {outputfile}...")
                try:
                    upload_file_to_zenodo(cachefile, job.rule.name, tarball=tarball)
                except exceptions.ZenodoException as e:
                    # NOTE: we treate all Zenodo caching errors as non-fatal
                    logger.error(str(e))

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
            try:
                upload_file_to_zenodo(cachefile, job.rule.name, tarball=tarball)
            except exceptions.ZenodoException as e:
                # NOTE: we treate all Zenodo caching errors as non-fatal
                logger.error(str(e))

        return result

    # Apply them
    output_file_cache.fetch = types.MethodType(fetch, output_file_cache)
    output_file_cache.store = types.MethodType(store, output_file_cache)


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
        if ur.conda_env:
            if hasattr(ur.conda_env, "is_file") and ur.conda_env.is_file:
                ur.set_input(str(ur.conda_env.file))
            else:
                ur.set_input(ur.conda_env)
        else:
            # All user rules should run in conda envs!
            # TODO
            raise exceptions.MissingCondaEnvironmentInUserRule()