from .logging import get_logger
import sys
from contextlib import contextmanager


@contextmanager
def no_traceback():
    """
    Sets a custom exception handler for the scope of a 'with' block.

    https://stackoverflow.com/a/40347369
    """
    sys.excepthook = no_traceback_excepthook
    yield
    sys.excepthook = sys.__excepthook__


def no_traceback_excepthook(type, value, traceback):
    """
    Print an exception message without the traceback.

    """
    print(": ".join([str(type.__name__), str(value)]))


class ShowyourworkException(Exception):
    def __init__(self, message="An error occurred while executing the workflow."):
        get_logger().error(message)
        super().__init__(message)


class MissingFigureOutputError(ShowyourworkException):
    pass


class FigureGenerationError(ShowyourworkException):
    pass


class TectonicError(ShowyourworkException):
    pass


class FigureFormatError(ShowyourworkException):
    pass


class MissingXMLFile(ShowyourworkException):
    pass


class GraphicsPathError(ShowyourworkException):
    pass


class ConfigError(ShowyourworkException):
    pass


class MissingConfigFile(ShowyourworkException):
    pass


class ZenodoError(ShowyourworkException):
    def __init__(self, data):
        data["status"] = data.get("status", "unknown")
        data["message"] = data.get(
            "message", "An error occurred while accessing Zenodo."
        )
        message = "Zenodo error {}: {}".format(data["status"], data["message"])
        super().__init__(message)


class MissingZenodoAccessToken(ShowyourworkException):
    def __init__(self, token_name):
        message = (
            f"Zenodo access token `{token_name}` not found. "
            "This should be set as an environment variable."
        )
        super().__init__(message)


class ZenodoRecordNotFound(ShowyourworkException):
    def __init__(self, record_id, id_type="version or concept"):
        message = (
            f"The provided `id` {record_id} does "
            "not seem to be a valid Zenodo {id_type} id."
        )
        super().__init__(message)


class ZenodoUploadError(ShowyourworkException):
    pass


class NotImplementedError(ShowyourworkException):
    pass


class TarballExtractionError(ShowyourworkException):
    pass