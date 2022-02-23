from .logging import get_logger
import sys
from contextlib import contextmanager


@contextmanager
def no_traceback():
    """
    Remove the traceback from an exception message.

    https://stackoverflow.com/a/63657211
    """
    tracebacklimit = getattr(sys, "tracebacklimit", 1000)
    sys.tracebacklimit = 0
    yield
    sys.tracebacklimit = tracebacklimit


class ShowyourworkException(Exception):
    def __init__(self, message="An error occurred while executing the workflow."):
        get_logger().error(message)
        super().__init__()


class ConfigError(ShowyourworkException):
    pass


class MissingFigureOutputError(ShowyourworkException):
    pass


class FigureGenerationError(ShowyourworkException):
    pass


class TectonicError(ShowyourworkException):
    def __init__(self, logfile=None):
        if logfile:
            with open(logfile, "r") as f:
                tectonic_log = f.readlines()

            # Ensure the user imported showyourwork
            for line in tectonic_log:
                if "Package: showyourwork" in line:
                    showyourwork_imported = True
                    break
            else:
                showyourwork_imported = False

            if showyourwork_imported:

                # Scan the log for an error message
                for i, line in enumerate(tectonic_log[::-1]):
                    if line.startswith("!"):
                        message = "".join(tectonic_log[-i - 1 :])
                        break
                else:
                    message = "An error occurred while compiling the manuscript."
                message += (
                    f"\nFor more information, check out the log file:\n{logfile}."
                )
            else:

                # Admonish the user (:
                message = r"Failed to compile manuscript. Did you forget to `\usepackage{showyourwork}`?"

        else:

            # No log to scan, so we're stuck with an uninformative message...
            message = "An error occurred while compiling the manuscript."

        super().__init__(message)


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
    def __init__(self, status="", message="An error occurred while accessing Zenodo."):
        super().__init__(f"Zenodo error {status}: {message}")


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


class ZenodoContentsError(ShowyourworkException):
    pass


class InvalidZenodoIdType(ShowyourworkException):
    pass


class TarballExtractionError(ShowyourworkException):
    pass


class MissingCondaEnvironmentInUserRule(ShowyourworkException):
    pass


class RunDirectiveNotAllowedInUserRules(ShowyourworkException):
    pass


class InvalidZenodoNotesField(ShowyourworkException):
    pass


class MissingOverleafCredentials(ShowyourworkException):
    pass


class MultipleOverleafIds(ShowyourworkException):
    def __init__(self):
        super().__init__(
            "Only a single Overleaf project ID may be specified in the config file."
        )


class OverleafError(ShowyourworkException):
    pass


class CalledProcessError(ShowyourworkException):
    pass


class OverleafAuthenticationError(ShowyourworkException):
    def __init__(self):
        super().__init__(
            "Overleaf authentication failed.\nMake sure you have set the environment "
            "variables and GitHub secrets `OVERLEAF_EMAIL` and `OVERLEAF_PASSWORD`. "
            "See the docs for details."
        )


# --


class FileNotFoundOnZenodo(Exception):
    """
    Note: not a subclass of `ShowyourworkException` since we
    don't want this printed to the logs.

    """

    def __init__(self, file_name):
        message = f"File {file_name} not found on Zenodo."
        super().__init__(message)