from . import exceptions
from .zenodo import download_file_from_zenodo, upload_file_to_zenodo
from .logging import get_logger
import time
import types
import os

try:
    import snakemake
except ModuleNotFoundError:
    snakemake = None


def patch_snakemake_cache(concept_id):
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
                    logger.info(
                        f"Searching Zenodo file cache: {outputfile}..."
                    )
                    download_file_from_zenodo(
                        cachefile,
                        job.rule.name,
                        concept_id,
                        tarball=tarball,
                    )
                    logger.info(
                        f"Restoring from Zenodo cache: {outputfile}..."
                    )
                except exceptions.FileNotFoundOnZenodo:
                    # Cache miss; not fatal
                    logger.warn(
                        f"Required version of file not found in cache: {outputfile}."
                    )
                    logger.warn(f"Running rule from scratch...")
                except exceptions.ZenodoException as e:
                    # NOTE: we treat all Zenodo caching errors as non-fatal
                    logger.error(str(e))
                    logger.warn(f"Running rule from scratch...")
            else:
                logger.info(f"Restoring from local cache: {outputfile}...")
                # Always ensure the cached file on Zenodo is the latest
                # local cache hit. (If it is, this is a no-op)
                logger.info(f"Syncing file with Zenodo cache: {outputfile}...")
                try:
                    upload_file_to_zenodo(
                        cachefile, job.rule.name, concept_id, tarball=tarball
                    )
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
                upload_file_to_zenodo(
                    cachefile, job.rule.name, concept_id, tarball=tarball
                )
            except exceptions.ZenodoException as e:
                # NOTE: we treate all Zenodo caching errors as non-fatal
                logger.error(str(e))

        return result

    # Apply them
    output_file_cache.fetch = types.MethodType(fetch, output_file_cache)
    output_file_cache.store = types.MethodType(store, output_file_cache)


def patch_snakemake_wait_for_files():
    """"""

    def wait_for_files(
        files, latency_wait=3, force_stay_on_remote=False, ignore_pipe=False
    ):
        """Wait for given files to be present in the filesystem."""
        files = list(files)

        def get_missing():
            return [
                f
                for f in files
                if not (
                    f.exists_remote
                    if (
                        isinstance(f, snakemake.io._IOFile)
                        and f.is_remote
                        and (force_stay_on_remote or f.should_stay_on_remote)
                    )
                    else os.path.exists(f)
                    if not (snakemake.io.is_flagged(f, "pipe") and ignore_pipe)
                    else True
                )
            ]

        missing = get_missing()
        if missing:
            get_logger().info(
                "Waiting at most {} seconds for missing files.".format(
                    latency_wait
                )
            )
            for _ in range(latency_wait):
                missing = get_missing()
                if not missing:
                    return
                time.sleep(1)
            missing = "\n".join(get_missing())
            raise exceptions.MissingFigureOutputError(
                f"Missing files after {latency_wait} seconds:\n" f"{missing}"
            )

    # Apply the patch
    snakemake.wait_for_files = wait_for_files
    snakemake.io.wait_for_files = wait_for_files
    snakemake.dag.wait_for_files = wait_for_files
    snakemake.jobs.wait_for_files = wait_for_files