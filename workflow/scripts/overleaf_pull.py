import sys


# Import utils
sys.path.insert(1, snakemake.config["workflow_abspath"])
from utils import overleaf


# Pull files from Overleaf
project_id = snakemake.config["overleaf"].get("id", None)
if project_id:
    files = snakemake.config["overleaf"].get("pull", [])
    if files:
        overleaf.pull_files(files, project_id)