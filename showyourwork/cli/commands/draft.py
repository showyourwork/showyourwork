from ... import paths
from .run_snakemake import run_snakemake


def draft(snakemake_args=[], cores=1, conda_frontend="conda"):
    """Build the article in draft mode.

    This function builds the article PDF by running ``Snakemake`` in an isolated
    conda environment in draft mode.

    Args:
        snakemake_args (list): Additional arguments to pass to ``Snakemake``.

    """
    snakefile = paths.showyourwork().workflow / "draft.smk"
    run_snakemake(
        snakefile.as_posix(),
        run_type="draft",
        cores=cores,
        conda_frontend=conda_frontend,
        extra_args=snakemake_args,
        check=True,
    )
