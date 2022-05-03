from showyourwork import exceptions
from showyourwork.zenodo import Zenodo, Sandbox
from showyourwork.logging import get_logger
import sys
import subprocess

if __name__ == "__main__":

    # Snakemake config (available automagically)
    config = snakemake.config  # type:ignore

    # Get params
    doi = snakemake.params["doi"]
    deposit_id = doi.split("/")[-1]
    remote_file = snakemake.params["remote_file"]
    output = snakemake.output[0]

    # Initialize the logger
    logger = get_logger()

    # Get service url
    if doi.startswith(Zenodo.doi_prefix):
        url = Zenodo.url
    else:
        url = Sandbox.url

    # Download it
    progress_bar = ["--progress-bar"] if not config["github_actions"] else []
    result = subprocess.run(
        [
            "curl",
            "-f",
            f"https://{url}/record/{deposit_id}/files/{remote_file}",
            *progress_bar,
            "--output",
            output,
        ]
    )
    if result.returncode != 0:
        raise exceptions.ZenodoDownloadError()