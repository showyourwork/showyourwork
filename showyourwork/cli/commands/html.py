from ... import paths
from .run_snakemake import run_snakemake


def html(snakemake_args=[], cores=1, conda_frontend="conda"):
    """Build the html.

    Args:
        options (str, optional): Additional options to pass to Snakemake.
    """
    snakefile = paths.showyourwork().workflow / "build.smk"
    run_snakemake(
        snakefile.as_posix(),
        run_type="html",
        cores=cores,
        conda_frontend=conda_frontend,
        extra_args=list(snakemake_args) + ["syw__ar5ivist_entrypoint"],
        check=True,
    )
