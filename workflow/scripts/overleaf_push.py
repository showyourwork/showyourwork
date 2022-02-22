import sys


# Import utils
sys.path.insert(1, snakemake.config["workflow_abspath"])
from utils import overleaf


# Push files to Overleaf
project_id = snakemake.config["overleaf"].get("id", None)
if project_id:
    files = snakemake.config["overleaf"].get("push", [])
    if files:
        overleaf.push_files(files, project_id)