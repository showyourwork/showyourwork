from .base import ShowyourworkException


class RequestError(ShowyourworkException):
    def __init__(
        self,
        status="",
        message="An error occurred while accessing a remote server.",
    ):
        super().__init__(f"Request error {status}: {message}")


class ConfigError(ShowyourworkException):
    pass


class MissingFigureOutputError(ShowyourworkException):
    pass


class FigureGenerationError(ShowyourworkException):
    pass


class ConfigError(ShowyourworkException):
    pass


class MissingConfigFile(ShowyourworkException):
    pass


class NotImplementedError(ShowyourworkException):
    pass


class TarballExtractionError(ShowyourworkException):
    pass


class MissingCondaEnvironmentInUserRule(ShowyourworkException):
    pass


class RunDirectiveNotAllowedInUserRules(ShowyourworkException):
    pass


class CalledProcessError(ShowyourworkException):
    pass
