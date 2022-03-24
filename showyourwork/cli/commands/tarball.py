from ... import paths
from ..conda_env import run
from ..error_handler import error_handler
from ..constants import SNAKEMAKE


def tarball(options=""):
    """Build the article tarball."""
    snakefile = paths.showyourwork().workflow / "build.smk"
    command = f"{SNAKEMAKE} {options} -s {snakefile} syw__arxiv_entrypoint"
    result = run(command, check=False)
    error_handler(result.returncode)