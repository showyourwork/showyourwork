import sys
from pathlib import Path
import subprocess


# Import utils
sys.path.insert(1, snakemake.config["workflow_abspath"])
from utils import exceptions, get_logger


# Get params
zenodo_url = snakemake.params["zenodo_url"]
deposit_id = snakemake.params["deposit_id"]
remote_file = snakemake.params["remote_file"]
output = snakemake.output[0]


# Initialize the logger
logger = get_logger()


# Download it
result = subprocess.run(
    [
        "curl",
        f"https://{zenodo_url}/record/{deposit_id}/files/{remote_file}",
        "--output",
        output,
    ]
)
if result.returncode != 0:
    # TODO
    raise exceptions.ZenodoError()

# TODO: Catch 404s