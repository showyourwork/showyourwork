from ... import paths
from .run_snakemake import run_snakemake


def tarball(snakemake_args=[], cores=1, conda_frontend="conda"):
    """Build the article tarball.

    Args:
        options (str, optional): Additional options to pass to Snakemake.
    """
    snakefile = paths.showyourwork().workflow / "build.smk"
    run_snakemake(
        snakefile.as_posix(),
        run_type="tarball",
        cores=cores,
        conda_frontend=conda_frontend,
        extra_args=list(snakemake_args) + ["syw__arxiv_entrypoint"],
        check=True,
    )
