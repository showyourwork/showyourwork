from pathlib import Path
import inspect
import logging
import snakemake


__all__ = ["setup_logging", "get_snakemake_variable"]


def setup_logging(verbose=False, logfile=None):

    # Get the default Snakemake logger
    logger = snakemake.logging.logger

    # Suppress all non-critical Snakemake output
    # to the terminal (unless verbose); save it all
    # for the logs!
    for handler in logger.logger.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.setLevel(logging.DEBUG)
        else:
            if not verbose:
                handler.setLevel(logging.ERROR)

    # Symlink the logfile in the showyourwork temp dir
    if logfile is not None:
        smlogfile = logger.get_logfile()
        smlogfile = Path(smlogfile)
        logfile = Path(logfile)
        if logfile.exists():
            logfile.unlink()
        logfile.symlink_to(smlogfile)


def get_snakemake_variable(name, default=None):
    """
    Infer the value of a variable within snakemake.

    This is extremely hacky.
    """
    for level in inspect.stack():
        value = level.frame.f_locals.get(name, None)
        if value is not None:
            return value
    return default