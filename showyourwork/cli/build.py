from .. import __version__
from .. import paths
from .conda_env import run
from .error_handler import error_handler


def _build(options=""):
    """Build the article."""
    # Preprocess, then main build
    precmd = f"SNAKEMAKE_OUTPUT_CACHE={paths.user().cache}"
    defaults = "-c1 --use-conda --reason --cache"
    for file in ["preprocess.smk", "build.smk"]:
        snakefile = paths.showyourwork().workflow / file
        command = f"{precmd} snakemake {defaults} {options} -s {snakefile}"
        result = run(command, check=False)
        error_handler(result.returncode)