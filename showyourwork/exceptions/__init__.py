"""
Defines custom exceptions for the ``showyourwork`` package.

"""

from .base import ShowyourworkException as ShowyourworkException
from .github import (
    GitHubException as GitHubException,
    MissingGitHubAPIKey as MissingGitHubAPIKey,
)
from .latex import (
    FigureFormatError as FigureFormatError,
    GraphicsPathError as GraphicsPathError,
    LaTeXException as LaTeXException,
    MissingXMLFile as MissingXMLFile,
    TectonicError as TectonicError,
    UnableToInferClassName as UnableToInferClassName,
)
from .other import (
    CondaNotFoundError as CondaNotFoundError,
    CondaVersionError as CondaVersionError,
    ConfigError as ConfigError,
    FigureGenerationError as FigureGenerationError,
    MissingConfigFile as MissingConfigFile,
    MissingDependencyError as MissingDependencyError,
    MissingFigureOutputError as MissingFigureOutputError,
    RequestError as RequestError,
    ShowyourworkNotFoundError as ShowyourworkNotFoundError,
)
from .overleaf import (
    MissingOverleafCredentials as MissingOverleafCredentials,
    MultipleOverleafIds as MultipleOverleafIds,
    OverleafAuthenticationError as OverleafAuthenticationError,
    OverleafError as OverleafError,
    OverleafException as OverleafException,
    OverleafRateLimitExceeded as OverleafRateLimitExceeded,
)
from .zenodo import (
    FileNotFoundOnZenodo as FileNotFoundOnZenodo,
    InvalidZenodoDOI as InvalidZenodoDOI,
    InvalidZenodoIdType as InvalidZenodoIdType,
    InvalidZenodoNotesField as InvalidZenodoNotesField,
    MissingZenodoAccessToken as MissingZenodoAccessToken,
    ZenodoContentsError as ZenodoContentsError,
    ZenodoDownloadError as ZenodoDownloadError,
    ZenodoError as ZenodoError,
    ZenodoException as ZenodoException,
    ZenodoRecordNotFound as ZenodoRecordNotFound,
    ZenodoUploadError as ZenodoUploadError,
)
