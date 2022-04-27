from .base import ShowyourworkException


class GitHubException(ShowyourworkException):
    pass


class MissingGitHubAPIKey(GitHubException):
    def __init__(self, token_name):
        message = (
            f"GitHub API key `{token_name}` not found. "
            "This should be set as both an environment variable "
            "and a GitHub repository secret."
        )
        super().__init__(message)