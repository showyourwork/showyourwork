__all__ = ["ShowyourworkException", "MissingOutputError"]


class ShowyourworkException(Exception):
    pass


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