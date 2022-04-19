from . import paths
from pathlib import Path
import logging
import platform
import threading
import os
import sys

try:
    import snakemake
except ModuleNotFoundError:
    snakemake = None


__all__ = ["get_logger", "setup_logging", "clear_errors"]


class ColorizingStreamHandler(logging.StreamHandler):
    """
    Adapted from snakemake.logging.ColorizingStreamHandler.

    Adds color support to stdout, if available.

    """

    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
    RESET_SEQ = "\033[0m"
    COLOR_SEQ = "\033[%dm"
    BOLD_SEQ = "\033[1m"
    colors = {
        "WARNING": YELLOW,
        "INFO": GREEN,
        "DEBUG": BLUE,
        "CRITICAL": RED,
        "ERROR": RED,
    }

    def __init__(self, stream=sys.stderr):
        super().__init__(stream=stream)
        self._output_lock = threading.Lock()
        self.nocolor = not self.can_color_tty()

    def can_color_tty(
        self,
    ):
        if "TERM" in os.environ and os.environ["TERM"] == "dumb":
            return False
        return self.is_tty and not platform.system() == "Windows"

    @property
    def is_tty(self):
        isatty = getattr(self.stream, "isatty", None)
        return isatty and isatty()

    def emit(self, record):
        with self._output_lock:
            try:
                self.format(record)
                self.stream.write(self.decorate(record))
                self.stream.write(getattr(self, "terminator", "\n"))
                self.flush()
            except BrokenPipeError as e:
                raise e
            except (KeyboardInterrupt, SystemExit):
                pass
            except Exception as e:
                self.handleError(record)

    def decorate(self, record):
        message = record.message
        message = [message]
        if not self.nocolor and record.levelname in self.colors:
            message.insert(
                0, self.COLOR_SEQ % (30 + self.colors[record.levelname])
            )
            message.append(self.RESET_SEQ)
        return "".join(message)


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
        return message


def get_logger():
    """
    Get our custom terminal logger.

    """
    logger = logging.getLogger("showyourwork")
    if not logger.handlers:

        # Root level
        logger.setLevel(logging.DEBUG)

        # Terminal: all messages
        stream_handler = ColorizingStreamHandler()
        stream_handler.setLevel(logging.INFO)
        logger.addHandler(stream_handler)

        try:

            LOGS = paths.user().logs

        except:

            # Can't resolve path to logs; assume we're not
            # in a showyourwork/git repo and fail silently.
            pass

        else:

            # File: just showyourwork errors
            error_file = LOGS / "showyourwork_errors.log"
            file_handler = logging.FileHandler(error_file)
            file_handler.setLevel(logging.ERROR)
            logger.addHandler(file_handler)

            # File: all showyourwork messages
            msg_file = LOGS / "showyourwork.log"
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
    if not snakemake:
        return

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
        snakemake_logger.custom_error_handler.setFormatter(
            SnakemakeFormatter()
        )
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