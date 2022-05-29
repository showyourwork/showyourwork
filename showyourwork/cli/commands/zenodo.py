from ...git import get_repo_branch
from ...zenodo import Zenodo
from ...config import edit_yaml
from ... import paths, exceptions, logging
import json
import os


def zenodo_publish(branch):
    """Publish the cache.

    Copies the latest draft on Zenodo Sandbox to Zenodo and publishes it,
    assigning it a permanent DOI.

    Args:
        branch (str): Branch to publish.
    """
    # Ensure there's a sandbox cache record for this branch
    if branch is None:
        branch = get_repo_branch()
    try:
        with edit_yaml("zenodo.yml") as config:
            sandbox_doi = config["cache"][branch]["sandbox"]
            assert sandbox_doi is not None
    except:
        raise exceptions.ShowyourworkException(
            f"Zenodo Sandbox deposit not found for branch {branch}."
        )
    sandbox = Zenodo(sandbox_doi)

    # Get the Zenodo doi, or create one if needed
    with edit_yaml("zenodo.yml") as config:
        zenodo_doi = config["cache"].get(branch, {}).get("zenodo", None)
        if "SANDBOX_ONLY" in os.environ:
            # Publish to Sandbox if we're in a test run
            service = "sandbox"
        else:
            service = "zenodo"
        if zenodo_doi is None:
            zenodo_doi = sandbox.copy_draft(service, branch=branch)
        else:
            zenodo_doi = sandbox.copy_draft(zenodo_doi, branch=branch)
    zenodo = Zenodo(zenodo_doi)

    # Publish the latest draft
    zenodo.publish()

    # Fill in the config
    with edit_yaml("zenodo.yml") as config:
        config["cache"][branch] = config["cache"].get(branch, {})
        config["cache"][branch]["zenodo"] = zenodo_doi


def zenodo_freeze(branch):
    """
    Freeze the Zenodo Sandbox deposit.

    This command publishes the latest draft on Zenodo Sandbox. This
    action makes the cache public and (semi-)permanent. Useful during
    development.

    Args:
        branch (str): Branch whose cache is to be frozen.
    """
    if branch is None:
        branch = get_repo_branch()
    try:
        with edit_yaml("zenodo.yml") as config:
            doi = config["cache"][branch]["sandbox"]
            assert doi is not None
    except:
        raise exceptions.ShowyourworkException(
            f"Zenodo Sandbox deposit not found for branch {branch}."
        )
    Zenodo(doi).publish()
    logging.get_logger().info(
        f"Zenodo Sandbox deposit draft with DOI {doi} published."
    )


def zenodo_create(branch):
    """
    Create a new Zenodo Sandbox deposit.

    This command creates a new draft on Zenodo Sandbox.

    Args:
        branch (str): Branch to create the cache for.
    """
    if branch is None:
        branch = get_repo_branch()
    try:
        with edit_yaml("zenodo.yml") as config:
            doi = config["cache"].get(branch, {}).get("sandbox", None)
            assert doi is None
    except AssertionError:
        raise exceptions.ShowyourworkException(
            f"Branch {branch} already has an associated Zenodo Sandbox deposit."
        )
    except:
        pass
    deposit = Zenodo("sandbox", branch=branch)
    with edit_yaml("zenodo.yml") as config:
        config["cache"] = config.get("cache", {})
        config["cache"][branch] = config["cache"].get(branch, {})
        config["cache"][branch]["sandbox"] = deposit.doi


def zenodo_delete(branch):
    """
    Delete the latest draft of the Zenodo Sandbox deposit.

    Args:
        branch (str): Branch whose cahce is to be deleted.
    """
    logger = logging.get_logger()
    if branch is None:
        branch = get_repo_branch()
    try:
        with edit_yaml("zenodo.yml") as config:
            doi = config["cache"][branch]["sandbox"]
            assert doi is not None
    except:
        raise exceptions.ShowyourworkException(
            f"Zenodo Sandbox deposit not found for branch {branch}."
        )
    Zenodo(doi).delete()
    with edit_yaml("zenodo.yml") as config:
        config["cache"][branch]["sandbox"] = None