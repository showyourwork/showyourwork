from . import paths
from pathlib import Path
import logging
import snakemake


__all__ = ["get_logger", "setup_logging", "clear_errors"]


def get_logger():
    """
    Get our custom terminal logger.

    """
    logger = logging.getLogger("showyourwork")
    if not logger.handlers:

        # Root level
        logger.setLevel(logging.DEBUG)

        # Terminal: all messages
        stream_handler = snakemake.logging.ColorizingStreamHandler()
        stream_handler.setLevel(logging.INFO)
        logger.addHandler(stream_handler)

        # File: just showyourwork errors
        error_file = paths.user().logs / "showyourwork_errors.log"
        file_handler = logging.FileHandler(error_file)
        file_handler.setLevel(logging.ERROR)
        logger.addHandler(file_handler)

        # File: all showyourwork messages
        msg_file = paths.user().logs / "showyourwork.log"
        file_handler = logging.FileHandler(msg_file)
        file_handler.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)

    return logger


def clear_snakemake_errors():
    error_file = paths.user().logs / "snakemake_errors.log"
    if error_file.exists():
        error_file.unlink()


def clear_showyourwork_errors():
    error_file = paths.user().logs / "showyourwork_errors.log"
    if error_file.exists():
        error_file.unlink()


def clear_errors():
    clear_snakemake_errors()
    clear_showyourwork_errors()

    # Remove temp checkpoint files
    for file in paths.user().checkpoints.glob("*"):
        file.unlink()


def setup_logging(verbose=False, logfile=None):
    """
    Hack the Snakemake logger to suppress most of its terminal output.

    """
    # Get our custom logger
    logger = get_logger()

    # Get the default Snakemake logger
    snakemake_logger = snakemake.logging.logger

    # Suppress *all* Snakemake output to the terminal (unless verbose);
    # save it all for the logs!
    for handler in snakemake_logger.logger.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.setLevel(logging.DEBUG)
        else:
            if not verbose:
                handler.setLevel(logging.CRITICAL)

    # Custom error handler
    if not hasattr(snakemake_logger, "custom_error_handler"):
        error_file = paths.user().logs / "snakemake_errors.log"
        snakemake_logger.custom_error_handler = logging.FileHandler(error_file)
        snakemake_logger.custom_error_handler.setLevel(logging.ERROR)
        snakemake_logger.logger.addHandler(
            snakemake_logger.custom_error_handler
        )

    # Log job messages to the terminal manually
    def job_info(self, **msg):
        msg["level"] = "job_info"
        self.handler(msg)
        if msg.get("msg", None):
            logger.info(msg["msg"])

    snakemake_logger.job_info = lambda **msg: job_info(snakemake_logger, **msg)

    # Allow all conda messages to come through
    snakemake.deployment.conda.logger = logger

    # Store the path to the snakemake logfile in the showyourwork temp dir
    # At the end of the build, we'll copy the contents over to this same file.
    if logfile is not None:
        smlogfile = snakemake_logger.get_logfile()
        if smlogfile:
            smlogfile = str(Path(smlogfile))
            with open(Path(logfile), "w") as f:
                print(smlogfile, file=f)
        else:
            # TODO: I've encountered cases where the logfile is None;
            # happened once when I had checkpoints in the workflow.
            # Investigate this, and whether we should change the logic here.
            pass