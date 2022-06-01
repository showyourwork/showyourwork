"""
Implements functions that modify the behavior of Snakemake.

This functionality is pretty hacky and is almost certainly not future-proof, 
so if we change the Snakemake version (currently pinned at 16.5.5), we may have 
to update the code in this file.

"""
from . import paths, exceptions
from .zenodo import Zenodo
from .logging import get_logger, ColorizingStreamHandler
import inspect
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


def patch_snakemake_missing_input_leniency():
    """
    Patches snakemake to raise an error if there are no producers for
    any of the output files in the DAG.

    """
    _dag_debug = snakemake.dag.logger.dag_debug

    def dag_debug(msg):
        if type(msg) is dict:
            if "No producers found, but file is present on disk" in msg.get(
                "msg", ""
            ):
                file = msg["file"]
                msg = str(msg["exception"])
                raise exceptions.MissingDependencyError(
                    f"File {file} is present on disk, but there "
                    f"is no valid rule to generate it.\n{msg}"
                )
        _dag_debug(msg)

    snakemake.dag.logger.dag_debug = dag_debug


def get_snakemake_variable(name, default=None):
    """
    Infer the value of a variable within snakemake.

    This is extremely hacky.
    """
    levels = inspect.stack()
    for level in levels:
        value = level.frame.f_locals.get(name, None)
        if value is not None:
            return value
    return default


def patch_snakemake_cache(zenodo_doi, sandbox_doi):
    """
    Patches the Snakemake cache functionality to

        - Add custom logging messages
        - Attempt to download the cache file from Zenodo or Zenodo Sandbox on `fetch()`
        - Uploads the cache file to Zenodo Sandbox on `store()`

    """
    # Get the showyourwork logger
    logger = get_logger()

    # The instance we'll patch
    output_file_cache = snakemake.workflow.workflow.output_file_cache

    if output_file_cache is not None:

        # Instantiate our interfaces
        if zenodo_doi is not None:
            zenodo = Zenodo(zenodo_doi)
        else:
            zenodo = None
        if sandbox_doi is not None:
            sandbox = Zenodo(sandbox_doi)
        else:
            sandbox = None

        # Make a copy of the original methods
        _fetch = output_file_cache.fetch
        _store = output_file_cache.store

        # Define the patches
        def fetch(self, job):
            # If the cache file is a directory, we must tar it up
            # Recall that cacheable jobs can only have a _single_ output
            # (unless using multiext), so checking the first output should suffice
            config = snakemake.workflow.config
            if job.output[0].is_directory:
                tarball = True
            else:
                tarball = False
            for outputfile, cachefile in self.get_outputfiles_and_cachefiles(
                job
            ):
                file_exists = cachefile.exists()
                if not file_exists:
                    # Attempt to download from Zenodo and then Zenodo Sandbox
                    try:
                        logger.info(
                            f"Searching remote file cache: {outputfile}..."
                        )
                        try:
                            logger.debug(f"Searching Zenodo cache...")
                            if zenodo is None:
                                logger.debug(
                                    f"Zenodo DOI not provided in `zenodo.yml`."
                                )
                                raise exceptions.FileNotFoundOnZenodo(None)
                            zenodo.download_file(
                                cachefile, job.rule.name, tarball=tarball
                            )
                            file_exists = True
                            logger.info(
                                f"Restoring from Zenodo cache: {outputfile}..."
                            )
                        except exceptions.FileNotFoundOnZenodo:
                            logger.debug(f"Searching Zenodo Sandbox cache...")
                            sandbox.download_file(
                                cachefile, job.rule.name, tarball=tarball
                            )
                            file_exists = True
                            logger.info(
                                f"Restoring from Zenodo Sandbox cache: {outputfile}..."
                            )
                    except exceptions.FileNotFoundOnZenodo:
                        # Cache miss; not fatal
                        exceptions.restore_trace()
                        logger.warn(
                            f"Required version of file not found in cache: {outputfile}."
                        )
                        if config.get("github_actions") and not config.get(
                            "run_cache_rules_on_ci"
                        ):
                            raise exceptions.ShowyourworkException(
                                f"Cache for {outputfile} not found, and "
                                "`run_cache_rules_on_ci` is set to `False`."
                            )
                        else:
                            logger.warn(f"Running rule from scratch...")
                    except exceptions.ZenodoException as e:
                        # NOTE: we treat all Zenodo caching errors as non-fatal
                        exceptions.restore_trace()
                        logger.error(str(e))
                        if config.get("github_actions") and not config.get(
                            "run_cache_rules_on_ci"
                        ):
                            raise exceptions.ShowyourworkException(
                                f"Cache for {outputfile} not found, and "
                                "`run_cache_rules_on_ci` is set to `False`."
                            )
                        else:
                            logger.warn(f"Running rule from scratch...")
                else:
                    # We're good
                    logger.info(f"Restoring from local cache: {outputfile}...")

                # If the file exists (either restored from the local cache or
                # from the Zenodo cache), we sync it back so that the latest
                # draft on Zenodo Sandbox is always fully up to date. Note that
                # we check the hash before uploading, so if it's already up to
                # date, this is a no-op. Note that GitHub Actions runs should
                # never update the cache
                if file_exists and not snakemake.workflow.config.get(
                    "github_actions"
                ):
                    logger.info(
                        f"Syncing file with Zenodo Sandbox cache: {outputfile}..."
                    )
                    try:
                        sandbox.upload_file(
                            cachefile,
                            job.rule.name,
                            tarball=tarball,
                        )
                    except exceptions.ZenodoException as e:
                        # NOTE: we treate all Zenodo caching errors as non-fatal
                        exceptions.restore_trace()
                        logger.error(str(e))

            # Call the original method
            return _fetch(job)

        def store(self, job):
            # Call the original method
            result = _store(job)

            # GitHub Actions runs should never update the cache
            if not snakemake.workflow.config.get("github_actions"):

                # See note in `fetch()` about tarballs
                if job.output[0].is_directory:
                    tarball = True
                else:
                    tarball = False

                for (
                    outputfile,
                    cachefile,
                ) in self.get_outputfiles_and_cachefiles(job):
                    logger.info(
                        f"Caching output file on remote: {outputfile}..."
                    )
                    try:
                        sandbox.upload_file(
                            cachefile, job.rule.name, tarball=tarball
                        )
                    except exceptions.ZenodoException as e:
                        # NOTE: we treate all Zenodo caching errors as non-fatal
                        exceptions.restore_trace()
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