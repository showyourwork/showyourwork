from . import paths, exceptions
from .zenodo import download_file_from_zenodo, upload_file_to_zenodo
from .logging import get_logger, ColorizingStreamHandler
import logging
import time
import types
import os

try:
    import snakemake
except ModuleNotFoundError:
    snakemake = None


class SnakemakeFormatter(logging.Formatter):
    """
    Format Snakemake errors before displaying them on stdout.

    Sometimes, Snakemake fails with suggestions for commands to fix certain
    issues. We intercept those suggestions here, replacing them with the
    corresponding showyourwork syntax for convenience.
    """

    replacements = {
        "snakemake --cleanup-metadata": "showyourwork build --cleanup-metadata",
        "rerun your command with the --rerun-incomplete flag": "rerun showyourwork build with the --rerun-incomplete flag",
        "It can be removed with the --unlock argument": "It can be removed by passing --unlock to showyourwork build",
    }

    def format(self, record):
        message = record.getMessage()
        for key, value in self.replacements.items():
            message = message.replace(key, value)
        record.message = message
        return message


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
                        zenodo_url="zenodo.org",
                    )
                    logger.info(
                        f"Restoring from Zenodo cache: {outputfile}..."
                    )
                except exceptions.FileNotFoundOnZenodo:
                    # Cache miss; not fatal
                    exceptions.enable_trace()
                    logger.warn(
                        f"Required version of file not found in cache: {outputfile}."
                    )
                    logger.warn(f"Running rule from scratch...")
                except exceptions.ZenodoException as e:
                    # NOTE: we treat all Zenodo caching errors as non-fatal
                    exceptions.enable_trace()
                    logger.error(str(e))
                    logger.warn(f"Running rule from scratch...")
            else:
                logger.info(f"Restoring from local cache: {outputfile}...")
                # Always ensure the cached file on Zenodo is the latest
                # local cache hit. (If it is, this is a no-op)
                logger.info(f"Syncing file with Zenodo cache: {outputfile}...")
                try:
                    upload_file_to_zenodo(
                        cachefile,
                        job.rule.name,
                        concept_id,
                        tarball=tarball,
                        zenodo_url="zenodo.org",
                    )
                except exceptions.ZenodoException as e:
                    # NOTE: we treate all Zenodo caching errors as non-fatal
                    exceptions.enable_trace()
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
                    cachefile,
                    job.rule.name,
                    concept_id,
                    tarball=tarball,
                    zenodo_url="zenodo.org",
                )
            except exceptions.ZenodoException as e:
                # NOTE: we treate all Zenodo caching errors as non-fatal
                exceptions.enable_trace()
                logger.error(str(e))

        return result

    # Apply them
    output_file_cache.fetch = types.MethodType(fetch, output_file_cache)
    output_file_cache.store = types.MethodType(store, output_file_cache)


def patch_snakemake_logging():
    """
    Hacks the Snakemake logger to suppress most of its terminal output,
    and redirects the rest to a custom log file.

    """
    if not snakemake:
        return

    # Get our custom logger
    logger = get_logger()

    # Customize the default Snakemake logger
    snakemake_logger = snakemake.logging.logger

    # Suppress *all* Snakemake output to the terminal (unless verbose);
    # save it all for the logs!
    for handler in snakemake_logger.logger.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.setLevel(logging.DEBUG)
        else:
            if not snakemake.workflow.config.get("verbose", False):
                handler.setLevel(logging.CRITICAL)

    # Custom Snakemake stdout handler
    if not hasattr(snakemake_logger, "custom_stream_handler"):
        snakemake_logger.custom_stream_handler = ColorizingStreamHandler()
        snakemake_logger.custom_stream_handler.setLevel(logging.ERROR)
        snakemake_logger.custom_stream_handler.setFormatter(
            SnakemakeFormatter()
        )
        snakemake_logger.logger.addHandler(
            snakemake_logger.custom_stream_handler
        )

    # Custom Snakemake file handler
    if not hasattr(snakemake_logger, "custom_file_handler"):
        try:

            LOGS = paths.user().logs

        except:

            # Can't resolve path to logs; assume we're not
            # in a showyourwork/git repo and fail silently.
            pass

        else:

            snakemake_logger.custom_file_handler = logging.FileHandler(
                paths.user().logs / "snakemake.log"
            )
            snakemake_logger.custom_file_handler.setLevel(logging.DEBUG)
            snakemake_logger.logger.addHandler(
                snakemake_logger.custom_file_handler
            )

    # Suppress Snakemake exceptions if we caught them on the
    # showyourwork side
    def job_error(self, **msg):
        msg["level"] = "job_error"
        try:
            FLAGS = paths.user().flags
        except:
            self.handler(msg)
        else:
            if (FLAGS / "DISABLE_SNAKEMAKE_EXCEPTIONS").exists():
                for handler in snakemake_logger.logger.handlers:
                    if not isinstance(handler, logging.FileHandler):
                        handler.setLevel(logging.CRITICAL)
            else:
                self.handler(msg)

    snakemake_logger.job_error = lambda **msg: job_error(
        snakemake_logger, **msg
    )

    # Allow job info messages to come through
    def job_info(self, **msg):
        msg["level"] = "job_info"
        self.handler(msg)
        if msg.get("msg", None):
            logger.info(msg["msg"])

    snakemake_logger.job_info = lambda **msg: job_info(snakemake_logger, **msg)

    # Allow all conda messages to come through
    snakemake.deployment.conda.logger = logger


def patch_snakemake_wait_for_files():
    """
    Replace Snakemake's `wait_for_files` method with a custom version.

    """

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