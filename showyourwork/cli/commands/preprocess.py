from ... import paths
from ..conda_env import run_in_env
from pathlib import Path
import os


def preprocess(snakemake_args=[]):
    """Pre-processing step for the article build."""
    snakefile = Path("${SYW_PATH}") / "workflow" / "prep.smk"
    snakemake = f"SNAKEMAKE_OUTPUT_CACHE={paths.user().cache} SNAKEMAKE_RUN_TYPE='preprocess' snakemake -c1 --use-conda --reason --cache"
    command = f"{snakemake} {' '.join(snakemake_args)} -s {snakefile}"
    result = run_in_env(command, check=False)
    if result.returncode > 0:
        os._exit(1)
