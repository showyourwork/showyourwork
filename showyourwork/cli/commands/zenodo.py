from ...git import get_repo_branch
from ...zenodo import Zenodo
from ...config import render_config
from ... import paths, exceptions, logging
import json


def zenodo_publish(doi):
    if doi is None:
        try:
            branch = get_repo_branch()
            config = render_config()
            doi = config["showyourwork"]["cache"][branch]
        except:
            raise exceptions.ShowyourworkException(
                "Unable to infer the current Zenodo deposit ID. "
                "Please provide it using the `--doi` option."
            )
    Zenodo(doi).publish()
    logging.get_logger().info(
        f"Zenodo deposit draft with DOI {doi} published."
    )


def zenodo_create(branch, service):
    if branch is None:
        branch = get_repo_branch()
    deposit = Zenodo(service, branch=branch)
    logging.get_logger().info(
        f"Please add the entry `{branch}: {deposit.doi}` to the `cache` section of the config file."
    )


def zenodo_delete(doi):
    logger = logging.get_logger()
    if doi is None:
        try:
            branch = get_repo_branch()
            config = render_config()
            doi = config["showyourwork"]["cache"][branch]
        except:
            raise exceptions.ShowyourworkException(
                "Unable to infer the current Zenodo deposit ID. "
                "Please provide it using the `--doi` option."
            )
    Zenodo(doi).delete()
    logging.get_logger().info(
        f"Zenodo deposit draft with DOI {doi} deleted. "
        "Please remove the corresponding `cache` entry from the config file."
    )