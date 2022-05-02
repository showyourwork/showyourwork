from ...git import get_repo_branch, get_repo_slug
from ...zenodo import delete_deposit, create_deposit, publish_deposit
from ... import paths, exceptions, logging
import json


def zenodo(publish_arg, create_arg, delete_arg):
    logger = logging.get_logger()
    if publish_arg != None:
        if publish_arg == "auto":
            try:
                with open(paths.user().temp / "config.json", "r") as f:
                    config = json.load(f)
                    concept_id = config["showyourwork"]["cache"]["zenodo"]
            except:
                raise exceptions.ShowyourworkException(
                    "Unable to infer the current Zenodo deposit ID. "
                    "Please provide it as an argument to `--delete`."
                )
        else:
            concept_id = publish_arg
        publish_deposit(concept_id)
        logger.info(
            f"Zenodo deposit {concept_id} published."
        )
    elif create_arg != None:
        if create_arg == "auto":
            branch = get_repo_branch()
        else:
            branch = create_arg
        slug = get_repo_slug()
        title = f"Data for {slug} [{branch}]"
        concept_id = create_deposit(title)
        logger.info(
            f"Zenodo deposit {concept_id} created. "
            "Please add this to the `cache` entry in the config file."
        )
    elif delete_arg != None:
        if delete_arg == "auto":
            try:
                with open(paths.user().temp / "config.json", "r") as f:
                    config = json.load(f)
                    concept_id = config["showyourwork"]["cache"]["zenodo"]
            except:
                raise exceptions.ShowyourworkException(
                    "Unable to infer the current Zenodo deposit ID. "
                    "Please provide it as an argument to `--delete`."
                )
        else:
            concept_id = delete_arg
        delete_deposit(concept_id)
        logger.info(
            f"Zenodo deposit {concept_id} deleted. "
            "Please remove the `cache` entry from the config file."
        )
