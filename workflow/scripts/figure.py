import subprocess
import sys
from snakemake.io import wait_for_files


# Import utils
sys.path.insert(1, snakemake.config["workflow_abspath"])
from utils import exceptions


# Run the command
result = subprocess.run(snakemake.params["command"], shell=True)
if result.returncode != 0:
    # TODO
    raise exceptions.FigureGenerationError()


# Preempt the snakemake "Missing files after X seconds"
# error by checking for them ourselves
try:
    wait_for_files(
        snakemake.output,
        latency_wait=snakemake.config["latency_wait"],
        force_stay_on_remote=not snakemake.config["assume_shared_fs"],
        ignore_pipe=True,
    )
except IOError as e:
    # TODO
    raise exceptions.MissingFigureOutputError()