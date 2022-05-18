from .base import ShowyourworkException


class OverleafException(ShowyourworkException):
    pass


class MultipleOverleafIds(OverleafException):
    def __init__(self):
        super().__init__(
            "Only a single Overleaf project ID may be specified in the config file."
        )


class OverleafError(OverleafException):
    pass


class OverleafRateLimitExceeded(OverleafException):
    pass


class MissingOverleafCredentials(OverleafException):
    def __init__(self, **kwargs):
        message = (
            f"Overleaf credentials `OVERLEAF_EMAIL` and/or "
            "`OVERLEAF_PASSWORD` not found. "
            "These should be set as both environment variables "
            "and GitHub repository secrets."
        )
        super().__init__(
            message,
            level=kwargs.get("level", "warn"),
        )


class OverleafAuthenticationError(OverleafException):
    def __init__(self, **kwargs):
        super().__init__(
            "Overleaf authentication failed. See the docs for details.",
            level=kwargs.get("level", "warn"),
        )
