"""
Downloads a publically available file from a Zenodo or Zenodo Sandbox record.

"""
from showyourwork import exceptions
from showyourwork.zenodo import Zenodo
from showyourwork.logging import get_logger
import sys
import subprocess

if __name__ == "__main__":

    # Snakemake config (available automagically)
    config = snakemake.config  # type:ignore

    # Get params
    doi = snakemake.params["doi"]
    remote_file = snakemake.params["remote_file"]
    output = snakemake.output[0]

    # Initialize the logger
    logger = get_logger()

    # Get the zenodo interface
    deposit = Zenodo(doi)

    # Download it
    progress_bar = ["--progress-bar"] if not config["github_actions"] else []
    result = subprocess.run(
        [
            "curl",
            "-f",
            f"https://{deposit.url}/record/{deposit.deposit_id}/files/{remote_file}",
            *progress_bar,
            "--output",
            output,
        ]
    )
    if result.returncode != 0:
        raise exceptions.ZenodoDownloadError()