import json
import os
from collections.abc import MutableMapping

import requests

from . import exceptions
from .subproc import parse_request


def flatten_dict(d, parent_key="", sep="__"):
    """
    Flatten a nested dictionary.

    Adapted from https://stackoverflow.com/a/6027615.
    """
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, MutableMapping):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


def get_access_token(token_name="GH_API_KEY", error_if_missing=False):
    """
    Return the access token stored in the environment variable `token_name`.

    Args:
        token_name (str, optional): Name of the environment variable containing
            the GitHub API token. Default is ``GH_API_KEY``.
        error_if_missing (bool, optional): Raise an error if the environment
            variable is not defined or empty? Default ``False``.

    """
    access_token = os.getenv(token_name, None)
    if error_if_missing and not access_token:
        raise exceptions.MissingGitHubAPIKey(token_name)
    return access_token


def get_authenticated_user():
    """
    Return the name of the current authenticated GitHub user.

    """
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
    Create a new repository on GitHub if it does not already exist.

    Args:
        name (str): The name of the repository (without the user).
        private (bool, optional): Set to ``True`` for private repos.
            Default ``False``.
        org (str, optional): Name of the organization (if applicable).
            Default is ``None``.

    """
    # Check if repo exists; if so, do nothing
    if org:
        url = f"https://api.github.com/repos/{org}/{name}"
    else:
        url = f"https://api.github.com/repos/{get_authenticated_user()}/{name}"
    if requests.get(url).status_code <= 204:
        return

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

    Args:
        name (str): The name of the repository (without the user).
        org (str, optional): Name of the organization (if applicable).
            Default is ``None``.
        quiet (bool, optional): Default ``False``.

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


def clear_cache(name, org=None):
    """
    Clear the Actions cache for a repository.

    Args:
        name (str): The name of the repository (without the user).
        org (str, optional): Name of the organization (if applicable).
            Default is ``None``.

    """
    if org:
        url = f"https://api.github.com/repos/{org}/{name}/actions/caches"
    else:
        url = f"https://api.github.com/repos/{get_authenticated_user()}/{name}/actions/caches"
    data = parse_request(
        requests.get(
            url,
            headers={
                "Accept": "application/vnd.github.v3+json",
                "Authorization": f"token {get_access_token()}",
            },
        )
    )
    for cache in data["actions_caches"]:
        parse_request(
            requests.delete(
                f"{url}/{cache['id']}",
                headers={
                    "Accept": "application/vnd.github.v3+json",
                    "Authorization": f"token {get_access_token()}",
                },
            )
        )


def get_workflow_run_status(name, org=None, q={}):
    """
    Checks the status of the latest GH Actions workflow run for a repository,
    optionally matching certain search criteria.

    Args:
        name (str): The name of the repository (without the user).
        org (str, optional): The name of the organization (if applicable).
            Default is ``None``.
        q (dict, optional): Search query, a dictionary containing all the
            key-value pairs to be matched in the JSON response of the
            workflow search. Nested dictionaries are allowed.

    Returns:
        tuple:
            A tuple containing ``(status, conclusion, url)``.

    The ``status`` is one of

    - ``queued``
    - ``in_progress``
    - ``completed``

    The ``conclusion`` is one of

    - ``action_required``
    - ``cancelled``
    - ``failure``
    - ``neutral``
    - ``success``
    - ``skipped``
    - ``stale``
    - ``timed_out``

    If the workflow is not found, both are set to ``unknown``.

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

    # Flatten the search criteria
    q = flatten_dict(q)

    # Return the first match
    for workflow_run in data["workflow_runs"]:

        # Flatten the response
        workflow_run = flatten_dict(workflow_run)

        # Check for matches
        for key, value in q.items():
            if workflow_run.get(key) != value:
                break
        else:
            return (
                workflow_run["status"],
                workflow_run["conclusion"],
                workflow_run["html_url"],
            )

    # No match
    return "unknown", "unknown", None
