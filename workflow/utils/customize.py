from . import paths
from pathlib import Path
import logging
import snakemake


def customize_logging(symlink=None):

    # Get the default Snakemake logger
    logger = snakemake.logging.logger

    # Suppress ALL non-critical Snakemake output
    # to the terminal; save it all for the logs!
    for handler in logger.logger.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.setLevel(logging.DEBUG)
        else:
            handler.setLevel(logging.CRITICAL)

    # Symlink the logfile in the showyourwork temp dir
    if symlink is not None:
        logfile = logger.get_logfile()
        symlink = Path(symlink)
        if symlink.exists():
            symlink.unlink()
        symlink.symlink_to(Path(logfile))