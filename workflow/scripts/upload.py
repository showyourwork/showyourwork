"""
This script is executed from within the ``zenodo`` rule
to upload files to Zenodo.

"""
import sys
from pathlib import Path

# Hack to import our custom functions
sys.path.insert(0, str(Path(__file__).parents[1]))
from helpers.exceptions import ShowyourworkException, ShowyourworkWarning
from helpers.zenodo import upload_simulation

# Go!
try:
    upload_simulation(
        snakemake.params["file_name"],
        snakemake.params["deposit_id"],
        snakemake.params["deposit_title"],
        snakemake.params["deposit_description"],
        snakemake.params["deposit_creators"],
        zenodo_url=snakemake.params["zenodo_url"],
        token_name=snakemake.params["token_name"],
        file_path=snakemake.params["file_path"],
        script=snakemake.params["script"],
        repo_url=snakemake.params["repo_url"],
    )
except ShowyourworkException as e:
    if "The server could not verify that you are authorized" in e.message:
        # Fail silently if the user doesn't have the
        # right authentication to upload to the deposit
        file_name = snakemake.params["file_name"]
        file_path = snakemake.params["file_path"]
        with open(f"{file_path}/{file_name}.zenodo", "w") as f:
            print("UNAUTHORIZED", file=f)
        # Raise a warning on job completion
        ShowyourworkWarning(
            e.message,
            script=e.script,
            rule_name=e.rule_name,
            context="This error occurred because showyourwork tried to "
            f"upload the file {file_name} to Zenodo under a deposit with id "
            f"{snakemake.params['deposit_id']}, but the current user does "
            "not have the correct authentication. The API token for Zenodo "
            f"should be stored in the environment variable {snakemake.params['token_name']}, "
            "but this token is either missing or invalid for editing the given "
            "deposit id. If you are a third-party user (i.e., you cloned someone "
            "else's repository and are attempting to build the paper locally), you "
            "can either run ``make fast`` "
            "to skip the dataset generation & upload step, or provide a different "
            "`id` to upload this file to in the config file; to reserve an `id` under "
            "your authentication, simply type `make reserve`. See the docs for more "
            "details.",
            brief=f"Unable to upload {file_name} to Zenodo.",
        )
    else:
        raise e