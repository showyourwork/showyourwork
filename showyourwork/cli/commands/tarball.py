from ... import paths
from ..conda_env import run_in_env
import os


def tarball(options=""):
    """Build the article tarball."""
    snakefile = paths.showyourwork().workflow / "build.smk"
    snakemake = f"SNAKEMAKE_OUTPUT_CACHE={paths.user().cache} snakemake -c1 --use-conda --reason --cache"
    command = f"{snakemake} {options} -s {snakefile} syw__arxiv_entrypoint"
    result = run_in_env(command, check=False)
    if result.returncode > 0:
        os._exit(1)
