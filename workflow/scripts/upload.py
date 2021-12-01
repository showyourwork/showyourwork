"""
This script is executed from within the ``zenodo`` rule
to upload files to Zenodo.

"""
import sys
from pathlib import Path

# Hack to import our custom functions
sys.path.insert(0, str(Path(__file__).parents[1]))
from helpers.exceptions import ShowyourworkException
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
        # Print the error, but don't make it fatal
        print(e.full_message)
    else:
        raise e