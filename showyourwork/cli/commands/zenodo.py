from ...git import get_repo_branch
from ...zenodo import Zenodo
from ...config import edit_yaml
from ... import paths, exceptions, logging
import json


def zenodo_publish(branch):
    if branch is None:
        branch = get_repo_branch()
    try:
        with edit_yaml("zenodo.yml") as config:
            doi = config["cache"][branch]
            assert doi is not None
    except:
        raise exceptions.ShowyourworkException(
            f"Zenodo deposit not found for branch {branch}."
        )
    Zenodo(doi).publish()
    logging.get_logger().info(
        f"Zenodo deposit draft with DOI {doi} published."
    )


def zenodo_create(branch, service):
    if branch is None:
        branch = get_repo_branch()
    try:
        with edit_yaml("zenodo.yml") as config:
            doi = config["cache"].get(branch, None)
            assert doi is None
    except AssertionError:
        raise exceptions.ShowyourworkException(
            f"Branch {branch} already has an associated Zenodo deposit."
        )
    except:
        pass
    deposit = Zenodo(service, branch=branch)
    with edit_yaml("zenodo.yml") as config:
        config["cache"][branch] = deposit.doi


def zenodo_delete(branch):
    logger = logging.get_logger()
    if branch is None:
        branch = get_repo_branch()
    try:
        with edit_yaml("zenodo.yml") as config:
            doi = config["cache"][branch]
            assert doi is not None
    except:
        raise exceptions.ShowyourworkException(
            f"Zenodo deposit not found for branch {branch}."
        )
    Zenodo(doi).delete()
    with edit_yaml("zenodo.yml") as config:
        config["cache"][branch] = None