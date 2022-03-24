from ... import paths
from ..conda_env import run
from ..error_handler import error_handler
from ..constants import SNAKEMAKE


def preprocess(options=""):
    """Pre-processing step for the article build."""
    snakefile = paths.showyourwork().workflow / "preprocess.smk"
    command = f"{SNAKEMAKE} {options} -s {snakefile}"
    result = run(command, check=False)
    error_handler(result.returncode)