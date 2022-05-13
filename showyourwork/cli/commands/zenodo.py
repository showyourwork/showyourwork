from ...git import get_repo_branch
from ...zenodo import Zenodo
from ...config import edit_yaml
from ... import paths, exceptions, logging
import json


def zenodo_publish(branch):
    if branch is None:
        branch = get_repo_branch()
    # TODO
    raise NotImplementedError("Coming soon!")


def zenodo_freeze(branch):
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
        config["cache"][branch] = config["cache"].get(branch, {})
        config["cache"][branch]["sandbox"] = deposit.doi


def zenodo_delete(branch):
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