from .base import ShowyourworkException


class RequestError(ShowyourworkException):
    def __init__(
        self,
        status="",
        message="An error occurred while accessing a remote server.",
    ):
        super().__init__(f"Request error {status}: {message}")


class CondaNotFoundError(ShowyourworkException):
    def __init__(self):
        super().__init__(
            f"Conda package manager not found. Is it installed and available in the system PATH?"
        )


class ShowyourworkNotFoundError(ShowyourworkException):
    def __init__(self, path):
        super().__init__(
            f"The requested version of showyourwork was not found at {path}."
        )


class ConfigError(ShowyourworkException):
    pass


class MissingFigureOutputError(ShowyourworkException):
    pass


class MissingDependencyError(ShowyourworkException):
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
    def __init__(self, name):
        super().__init__(
            f"The `run` directive is not allowed in user-defined rules. "
            f"Please use `script` or `shell` instead in rule {name}."
        )


class CalledProcessError(ShowyourworkException):
    pass
