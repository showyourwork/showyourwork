from .. import paths
from ..logging import get_logger
import traceback
import sys


def redirect_exception(*args, **kwargs):
    """
    Redirect the traceback exception printout to the log file.

    """
    exc = traceback.format_exception(*args, **kwargs)
    exc = "".join(exc)
    get_logger().debug(exc)


def custom_excepthook(cls, exc, tb):
    """
    Redirect the exception to the log file.

    """
    get_logger().debug("".join(traceback.format_exception(cls, exc, tb)))


print_exception = traceback.print_exception
excepthook = sys.excepthook


def disable_trace():
    """Disable the traceback from being printed to the screen.

    The traceback gets logged to file, unless the logging level
    is `DEBUG`, in which case it also gets printed to the screen.
    """
    traceback.print_exception = redirect_exception
    sys.excepthook = custom_excepthook
    try:
        flag = paths.user().flags / "DISABLE_SNAKEMAKE_EXCEPTIONS"
    except:
        pass
    else:
        flag.touch()


def restore_trace():
    """Restore traceback printing to the screen."""
    traceback.print_exception = print_exception
    sys.excepthook = excepthook
    try:
        flag = paths.user().flags / "DISABLE_SNAKEMAKE_EXCEPTIONS"
    except:
        pass
    else:
        flag.unlink()


class ShowyourworkException(Exception):
    def __init__(
        self,
        message="An error occurred while executing the workflow.",
        level="error",
    ):
        # Print the message using the logger
        if level == "error":
            get_logger().error(message)
        elif level == "warn":
            get_logger().warn(message)
        elif level == "info":
            get_logger().info(message)
        elif level == "debug":
            get_logger().debug(message)
        else:
            super().__init__(message)

        # Disable tracebacks; if this exception is caught, the
        # caller MUST call `restore_trace` in the `except` block
        # to restore traceback printing!
        disable_trace()

        # Raise the exception
        super().__init__()
