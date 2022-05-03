from showyourwork import exceptions
from showyourwork.zenodo import services
from showyourwork.logging import get_logger
import sys
import subprocess

if __name__ == "__main__":

    # Snakemake config (available automagically)
    config = snakemake.config  # type:ignore

    # Get params
    doi = snakemake.params["doi"]
    doi_prefix, deposit_id = doi.split("/")
    remote_file = snakemake.params["remote_file"]
    output = snakemake.output[0]

    # Initialize the logger
    logger = get_logger()

    # Get service url
    for service in services.values():
        if doi_prefix == service["doi_prefix"]:
            url = service["url"]
            break
    else:
        raise exceptions.ZenodoDownloadError()

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