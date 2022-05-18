from ... import paths
from ..conda_env import run_in_env
from pathlib import Path
import os


def tarball(options=""):
    """Build the article tarball."""
    snakefile = snakefile = Path("${SYW_PATH}") / "workflow" / "build.smk"
    snakemake = f"SNAKEMAKE_OUTPUT_CACHE={paths.user().cache} SNAKEMAKE_RUN_TYPE='tarball' snakemake -c1 --use-conda --reason --cache"
    command = f"{snakemake} {options} -s {snakefile} syw__arxiv_entrypoint"
    result = run_in_env(command, check=False)
    if result.returncode > 0:
        os._exit(1)
