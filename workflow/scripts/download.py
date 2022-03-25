from showyourwork import exceptions
from showyourwork.logging import get_logger
import sys
import subprocess


# Snakemake config (available automagically)
config = snakemake.config  # type:ignore


# Get params
zenodo_url = snakemake.params["zenodo_url"]
deposit_id = snakemake.params["deposit_id"]
remote_file = snakemake.params["remote_file"]
output = snakemake.output[0]


# Initialize the logger
logger = get_logger()


# Download it
progress_bar = ["--progress-bar"] if not config["github_actions"] else []
result = subprocess.run(
    [
        "curl",
        f"https://{zenodo_url}/record/{deposit_id}/files/{remote_file}",
        *progress_bar,
        "--output",
        output,
    ]
)
if result.returncode != 0:
    # TODO
    raise exceptions.ZenodoError()

# TODO: Catch 404s