from showyourwork.testing import run_snakemake


def test_config_args() -> None:
    run_snakemake(
        "tests/projects/plugins/staging/config",
        "staging2/stage.upload",
        "-s",
        "Snakefile_args",
        "--config",
        "working_directory=staging2",
    )


def test_config_kwargs() -> None:
    run_snakemake(
        "tests/projects/plugins/staging/config",
        "staging2/stage.upload",
        "-s",
        "Snakefile_kwargs",
        "--config",
        "working_directory=staging2",
    )
