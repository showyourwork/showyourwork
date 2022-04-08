from ..logging import get_logger
import sys
from contextlib import contextmanager


@contextmanager
def no_traceback():
    """
    Remove the traceback from an exception message.

    See https://stackoverflow.com/a/63657211
    """
    tracebacklimit = getattr(sys, "tracebacklimit", 1000)
    sys.tracebacklimit = 0
    yield
    sys.tracebacklimit = tracebacklimit


class ShowyourworkException(Exception):
    def __init__(
        self,
        message="An error occurred while executing the workflow.",
        level="error",
    ):
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
        super().__init__()