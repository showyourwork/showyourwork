from showyourwork.testing import run_snakemake


def test_noop_snapshot() -> None:
    run_snakemake("tests/projects/plugins/staging/noop-snapshot", "staging__upload")


def test_noop_restore() -> None:
    run_snakemake(
        "tests/projects/plugins/staging/noop-restore",
        "output/a.txt",
        "--config",
        "restore=True",
    )
