import os
import subprocess

from ... import paths


def build(snakemake_args=[]):
    """Build the article.

    This function builds the article PDF by running ``Snakemake`` in an isolated
    conda environment.

    Args:
        snakemake_args (list): Additional arguments to pass to ``Snakemake``.

    """
    snakefile = paths.showyourwork().workflow / "build.smk"
    snakemake = f"SNAKEMAKE_OUTPUT_CACHE={paths.user().cache} SNAKEMAKE_RUN_TYPE='build' snakemake -c1 --use-conda --conda-frontend conda --reason --cache"
    command = f"{snakemake} {' '.join(snakemake_args)} -s {snakefile}"
    result = subprocess.run(command, shell=True, check=False)
    if result.returncode > 0:
        os._exit(1)
