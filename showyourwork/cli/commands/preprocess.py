import os
from pathlib import Path

from ... import paths
from ..conda_env import run_in_env
from ..patches import SNAKEMAKE


def preprocess(snakemake_args=[]):
    """Pre-processing step for the article build.

    Args:
        snakemake_args (list, optional): Additional options to pass to Snakemake.
    """
    snakefile = Path("${SYW_PATH}") / "workflow" / "prep.smk"
    command_pre = f"SNAKEMAKE_OUTPUT_CACHE={paths.user().cache} SNAKEMAKE_RUN_TYPE='preprocess' {SNAKEMAKE} -c1 --use-conda --reason --cache"
    command = f"{command_pre} {' '.join(snakemake_args)} -s {snakefile}"
    result = run_in_env(command, check=False)
    if result.returncode > 0:
        os._exit(1)
