from snakemake.cli import parse_args

from ... import paths
from .run_snakemake import run_snakemake


def preprocess(snakemake_args=(), cores=1, conda_frontend="conda"):
    """Pre-processing step for the article build.

    Args:
        snakemake_args (list, optional): Additional options to pass to Snakemake.
    """
    snakemake_args = tuple(
        arg for arg in snakemake_args if arg not in ["--dry-run", "-n"]
    )
    snakefile = paths.showyourwork().workflow / "prep.smk"
    target_args = parse_args(snakemake_args)[1].targets
    snakemake_args_notargets = []
    for arg in snakemake_args:
        if arg not in target_args:
            snakemake_args_notargets.append(arg)
    run_snakemake(
        snakefile.as_posix(),
        run_type="preprocess",
        cores=cores,
        conda_frontend=conda_frontend,
        extra_args=snakemake_args_notargets,
        check=True,
    )
