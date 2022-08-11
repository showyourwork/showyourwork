import logging
import os
import platform
import sys
import threading

from . import paths

try:
    import snakemake
except ModuleNotFoundError:
    snakemake = None


__all__ = ["get_logger"]


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


def get_logger():
    """
    Return the custom showyourwork logger.

    Sets up the logging if needed.

    """
    logger = logging.getLogger("showyourwork")

    # Add showyourwork stream & file handlers
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

            # File: all showyourwork messages
            msg_file = LOGS / "showyourwork.log"
            file_handler = logging.FileHandler(msg_file)
            file_handler.setLevel(logging.DEBUG)
            logger.addHandler(file_handler)

    return logger
