import subprocess
import sys
from pathlib import Path
import time


# Import utils
sys.path.insert(1, snakemake.config["workflow_abspath"])
from utils import exceptions, get_logger, overleaf


# Initialize the logger
logger = get_logger()


# Run the command
result = subprocess.run(snakemake.params["command"], shell=True)
if result.returncode != 0:
    # TODO
    raise exceptions.FigureGenerationError()


# Preempt the snakemake "Missing files after X seconds"
# error by checking for them ourselves
latency_wait = snakemake.config["latency_wait"]
get_missing = lambda: [f for f in snakemake.output if not Path(f).exists()]
missing = get_missing()
if missing:
    logger.info("Waiting at most {} seconds for missing files.".format(latency_wait))
    for _ in range(latency_wait):
        if not get_missing():
            break
        time.sleep(1)
    # TODO
    raise exceptions.MissingFigureOutputError()


# Push new figures to Overleaf
project_id = snakemake.config["overleaf"].get("id", None)
if project_id:
    overleaf.push_files(snakemake.output, project_id)
