from ... import paths
from ..conda_env import run_in_env
from pathlib import Path
import os


def build(snakemake_args=[]):
    """Build the article.

    This function builds the article PDF by running ``Snakemake`` in an isolated
    conda environment.

    Args:
        snakemake_args (list): Additional arguments to pass to ``Snakemake``.

    """
    snakefile = Path("${SYW_PATH}") / "workflow" / "build.smk"
    snakemake = f"SNAKEMAKE_OUTPUT_CACHE={paths.user().cache} SNAKEMAKE_RUN_TYPE='build' snakemake -c1 --use-conda --reason --cache"
    command = f"{snakemake} {' '.join(snakemake_args)} -s {snakefile}"
    result = run_in_env(command, check=False)
    if result.returncode > 0:
        os._exit(1)
