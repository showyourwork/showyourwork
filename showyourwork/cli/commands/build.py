from ... import paths
from .run_snakemake import run_snakemake


def build(snakemake_args=[], cores=1, conda_frontend="conda"):
    """Build the article.

    This function builds the article PDF by running ``Snakemake`` in an isolated
    conda environment.

    Args:
        snakemake_args (list): Additional arguments to pass to ``Snakemake``.

    """
    snakefile = paths.showyourwork().workflow / "build.smk"
    run_snakemake(
        snakefile.as_posix(),
        run_type="build",
        cores=cores,
        conda_frontend=conda_frontend,
        extra_args=snakemake_args,
        check=True,
    )
