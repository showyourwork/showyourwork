from ... import paths
from ..conda_env import run
from ..error_handler import error_handler
from ..constants import SNAKEMAKE


def build(options=""):
    """Build the article."""
    snakefile = paths.showyourwork().workflow / "build.smk"
    command = f"{SNAKEMAKE} {options} -s {snakefile}"
    result = run(command, check=False)
    error_handler(result.returncode)