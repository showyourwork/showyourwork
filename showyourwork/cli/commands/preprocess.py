from ... import paths
from ..conda_env import run_in_env
from ..error_handler import error_handler


def preprocess(snakemake_args=[]):
    """Pre-processing step for the article build."""
    snakefile = paths.showyourwork().workflow / "preprocess.smk"
    snakemake = f"SNAKEMAKE_OUTPUT_CACHE={paths.user().cache} snakemake -c1 --use-conda --reason --cache"
    command = f"{snakemake} {' '.join(snakemake_args)} -s {snakefile}"
    result = run_in_env(command, check=False)
    error_handler(result.returncode)