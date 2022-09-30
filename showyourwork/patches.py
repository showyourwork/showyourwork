"""
Implements functions that modify the behavior of Snakemake.

This functionality is pretty hacky and is almost certainly not future-proof,
so if we change the Snakemake version (currently pinned at 16.5.5), we may have
to update the code in this file.

"""
import logging
import types

from . import exceptions, paths
from .logging import ColorizingStreamHandler, get_logger
from .zenodo import Zenodo

try:
    import snakemake
except ModuleNotFoundError:
    snakemake = None


def patch_snakemake_cache(zenodo_doi, sandbox_doi):
    """
    Patches the Snakemake cache functionality to

    - Add custom logging messages
    - Attempt to download the cache file from Zenodo/Zenodo Sandbox on ``fetch()``
    - Uploads the cache file to Zenodo Sandbox on ``store()``

    Args:
        zenodo_doi (str): The Zenodo DOI for the cache. Can be ``None``.
        sandbox_doi (str): The Zenodo Sandbox doi for the development cache.
            Can be ``None``.

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
                    except Exception as e:
                        # NOTE: we treat all Zenodo caching errors as non-fatal
                        exceptions.restore_trace()
                        logger.warn(
                            "File not found on remote cache. See logs for details."
                        )
                        if len(str(e)):
                            logger.debug(str(e))
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
                    except Exception as e:
                        # NOTE: we treate all Zenodo caching errors as non-fatal
                        exceptions.restore_trace()
                        logger.warn(
                            f"Failed to sync {outputfile} with Zenodo Sandbox cache. See logs for details."
                        )
                        if len(str(e)):
                            logger.debug(str(e))

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
                    except Exception as e:
                        # NOTE: we treate all Zenodo caching errors as non-fatal
                        exceptions.restore_trace()
                        logger.warn(
                            f"Failed to upload {outputfile} to Zenodo Sandbox cache. See logs for details."
                        )
                        if len(str(e)):
                            logger.debug(str(e))

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


def job_is_cached(job):
    """
    Return True if a job's outputs will be restored from cache.

    """
    logger = get_logger()
    cache = snakemake.workflow.workflow.output_file_cache

    # Check if user requested caching for job
    if not snakemake.workflow.workflow.is_cached_rule(job.rule):
        logger.debug(f"Job {job.name} is not cacheable.")
        return False

    # Check for a local cache hit
    try:
        if cache.exists(job):
            logger.debug(f"Local cache exists for job {job.name}.")
            return True
        else:
            pass
    except:
        # Job is not cacheable (no output files or multiple output files)
        logger.debug(f"Job {job.name} is not cacheable.")
        return False

    # Check if a remote cache exists
    branch = snakemake.workflow.config["git_branch"]
    zenodo_doi = snakemake.workflow.config["cache"][branch]["zenodo"]
    sandbox_doi = snakemake.workflow.config["cache"][branch]["sandbox"]
    if not zenodo_doi and not sandbox_doi:
        logger.debug(f"No cache hits for job {job.name}.")
        return False

    # Loop over cache files for the job (should really only be one)
    for outputfile, cachefile in cache.get_outputfiles_and_cachefiles(job):

        def file_exists(doi):
            """
            Return True if the file exists in the deposit with given doi.

            """
            if doi:
                try:
                    Zenodo(doi).download_file(
                        cachefile, job.rule.name, dry_run=True
                    )
                    return True
                except exceptions.FileNotFoundOnZenodo:
                    pass
            return False

        if file_exists(zenodo_doi) or file_exists(sandbox_doi):
            continue
        else:
            # No hit!
            logger.debug(f"No cache hits for job {job.name}.")
            return False

    # All cache files exist
    logger.debug(f"Remote cache exists for job {job.name}.")
    return True


def get_skippable_jobs(dag):
    """
    Search the DAG and return jobs we can safely skip due to
    downstream cache hits.

    """
    logger = get_logger()

    # Get all jobs (nodes) with cache hits
    cached_jobs = set([job for job in dag.jobs if job_is_cached(job)])
    logger.debug("The following jobs have cache hits:")
    logger.debug("    " + " ".join([job.name for job in cached_jobs]))

    # Nodes we can skip
    nodes = set(cached_jobs)

    # Loop until there are no changes to the graph
    new_nodes = True
    while new_nodes:

        new_nodes = set()
        for node in nodes:

            # Find all parents that are in `nodes`
            parents = set()
            for file in node.input:
                try:
                    parents |= set(dag.file2jobs(file)).difference(nodes)
                except snakemake.exceptions.MissingRuleException:
                    # Not all files have producers
                    pass

            # Find parents whose children are all in `nodes`
            for parent in parents:
                children = set()
                for file in parent.output:
                    children |= set(
                        [job for job in dag.jobs if file in job.input]
                    )
                if all([child in nodes for child in children]):
                    new_nodes.add(parent)

        # Add the new batch to the list of all skippable nodes
        nodes |= new_nodes

    # Exclude the cached jobs from the list we'll skip
    nodes = nodes.difference(cached_jobs)

    return nodes


def patch_snakemake_cache_optimization(dag):
    """
    Remove unnecessary jobs upstream of those with cache hits.

    See the full discussion about this feature
    `here <https://github.com/showyourwork/showyourwork/issues/124>`__.

    """
    logger = get_logger()

    if snakemake.workflow.workflow.output_file_cache is not None:

        # Get list of jobs we can skip
        skippable_jobs = get_skippable_jobs(dag)
        logger.debug(
            "Skipping the following jobs because of downstream cache hits:"
        )
        logger.debug("    " + " ".join([job.name for job in skippable_jobs]))

        # Patch the Snakemake method that executes jobs
        scheduler = snakemake.workflow.workflow.scheduler
        for executor in scheduler._executor, scheduler._local_executor:
            if executor:

                # Original method
                _cached_or_run = executor.cached_or_run

                # Intercept jobs we don't need to run
                def wrapper(self, job, run_func, *args):

                    if job in skippable_jobs:
                        # Instead of running this job, we create empty temp
                        # outputs to trick Snakemake & keep going
                        for output in job.output:
                            if not output.exists:
                                logger.warn(
                                    f"Skipping job {job.name} because "
                                    "of a downstream cache hit."
                                )
                                output.set_flags({"temp": True})
                                output.touch_or_create()
                        return
                    else:
                        # We need to run this rule
                        _cached_or_run(job, run_func, *args)

                # Patch the method
                executor.cached_or_run = types.MethodType(wrapper, executor)
