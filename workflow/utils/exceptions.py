from .logging import get_logger


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