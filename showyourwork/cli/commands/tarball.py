import os
from pathlib import Path

from ... import paths
from ..conda_env import run_in_env
from ..patches import SNAKEMAKE


def tarball(options=""):
    """Build the article tarball.

    Args:
        options (str, optional): Additional options to pass to Snakemake.
    """
    snakefile = snakefile = Path("${SYW_PATH}") / "workflow" / "build.smk"
    command_pre = f"SNAKEMAKE_OUTPUT_CACHE={paths.user().cache} SNAKEMAKE_RUN_TYPE='tarball' {SNAKEMAKE} -c1 --use-conda --reason --cache"
    command = f"{command_pre} {options} -s {snakefile} syw__arxiv_entrypoint"
    result = run_in_env(command, check=False)
    if result.returncode > 0:
        os._exit(1)
