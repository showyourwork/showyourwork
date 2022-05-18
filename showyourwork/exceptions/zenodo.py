from .base import ShowyourworkException


class ZenodoException(ShowyourworkException):
    pass


class ZenodoError(ZenodoException):
    def __init__(
        self, status="", message="An error occurred while accessing Zenodo."
    ):
        super().__init__(
            f"Zenodo error{' ' * bool(len(str(status)))}{status}: {message}"
        )


class ZenodoDownloadError(ZenodoException):
    def __init__(self):
        super().__init__(
            "An error occurred while downloading the dataset from Zenodo."
        )


class MissingZenodoAccessToken(ZenodoException):
    def __init__(self, token_name):
        message = (
            f"Zenodo access token `{token_name}` not found. "
            "This should be set as both an environment variable "
            "and a GitHub repository secret."
        )
        super().__init__(message)


class ZenodoRecordNotFound(ZenodoException):
    def __init__(self, record_id, id_type="version or concept"):
        message = (
            f"The provided `id` {record_id} does "
            f"not seem to be a valid Zenodo {id_type} id."
        )
        super().__init__(message)


class InvalidZenodoDOI(ZenodoException):
    def __init__(self, doi):
        message = (
            f"The provided `doi` {doi} does "
            f"not seem to be a valid Zenodo DOI."
        )
        super().__init__(message)


class ZenodoUploadError(ZenodoException):
    def __init__(self):
        super().__init__(
            "An error occurred while uploading the dataset from Zenodo."
        )


class ZenodoContentsError(ZenodoException):
    pass


class InvalidZenodoIdType(ZenodoException):
    pass


class InvalidZenodoNotesField(ZenodoException):
    def __init__(self):
        super().__init__(
            "The `Additional Notes` field in the current Zenodo record is not valid JSON."
        )


class FileNotFoundOnZenodo(ZenodoException):
    def __init__(self, file_name):
        super().__init__(f"File {file_name} not found on Zenodo.", level=None)