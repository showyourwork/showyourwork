from .subproc import parse_request
from . import exceptions
import requests
import os
import json


def get_access_token(token_name="GH_API_KEY", error_if_missing=False):
    """
    Return the access token stored in the environment variable `token_name`.

    """
    access_token = os.getenv(token_name, None)
    if error_if_missing and not access_token:
        raise exceptions.MissingGitHubAPIKey(token_name)
    return access_token


def get_authenticated_user():
    data = parse_request(
        requests.get(
            "https://api.github.com/user",
            headers={
                "Accept": "application/vnd.github.v3+json",
                "Authorization": f"token {get_access_token()}",
            },
        )
    )
    return data["login"]


def create_repo(name, description=None, private=False, org=None):
    """
    Create a new repository on GitHub.

    """
    # Delete repo (if it exists)
    delete_repo(name, org=org, quiet=True)

    # Create a new repo
    if description is None:
        description = "Nothing to see here!"
    if org:
        url = f"https://api.github.com/orgs/{org}/repos"
    else:
        url = "https://api.github.com/user/repos"
    parse_request(
        requests.post(
            url,
            data=json.dumps(
                {
                    "name": name,
                    "description": description,
                    "has_issues": False,
                    "has_projects": False,
                    "has_wiki": False,
                    "auto_init": False,
                    "private": private,
                }
            ),
            headers={
                "Accept": "application/vnd.github.v3+json",
                "Authorization": f"token {get_access_token()}",
            },
        )
    )


def delete_repo(name, org=None, quiet=False):
    """
    Delete a repository on GitHub.

    """
    if org:
        url = f"https://api.github.com/repos/{org}/{name}"
    else:
        url = f"https://api.github.com/repos/{get_authenticated_user()}/{name}"
    result = requests.delete(
        url,
        headers={
            "Accept": "application/vnd.github.v3+json",
            "Authorization": f"token {get_access_token()}",
        },
    )
    if not quiet:
        parse_request(result)


def get_latest_workflow_run_status(name, org=None):
    """
    Returns the tuple (status, conclusion, url).

    The `status` is one of

        queued, in_progress, completed

    The `conclusion` is one of

        action_required, cancelled, failure, neutral,
        success, skipped, stale, timed_out

    If the workflow is not found, both are set to `unknown`.

    The url is that of the workflow on GitHub actions.
    """
    if org:
        url = f"https://api.github.com/repos/{org}/{name}/actions/runs"
    else:
        url = f"https://api.github.com/repos/{get_authenticated_user()}/{name}/actions/runs"
    data = parse_request(
        requests.get(
            url,
            headers={
                "Accept": "application/vnd.github.v3+json",
                "Authorization": f"token {get_access_token()}",
            },
        )
    )
    if data["total_count"] == 0:
        return "unknown", "unknown", None

    # Get the first hit. _Should_ be the most recent workflow run, but
    # I don't know if that's guaranteed.
    workflow_run = data["workflow_runs"][0]
    return (
        workflow_run["status"],
        workflow_run["conclusion"],
        workflow_run["html_url"],
    )
