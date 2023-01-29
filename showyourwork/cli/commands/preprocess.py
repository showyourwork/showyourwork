from ... import paths
from .run_snakemake import run_snakemake


def preprocess(snakemake_args=[], cores=1, conda_frontend="conda"):
    """Pre-processing step for the article build.

    Args:
        snakemake_args (list, optional): Additional options to pass to Snakemake.
    """
    snakefile = paths.showyourwork().workflow / "prep.smk"
    run_snakemake(
        snakefile.as_posix(),
        run_type="preprocess",
        cores=cores,
        conda_frontend=conda_frontend,
        extra_args=snakemake_args,
        check=True,
    )
