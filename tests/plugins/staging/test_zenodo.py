import os
from typing import Generator

import pytest
from showyourwork.plugins.staging.zenodo_mock import ZenodoMock
from showyourwork.testing import run_snakemake


@pytest.fixture(scope="session")
def server() -> Generator[ZenodoMock, None, None]:
    server = ZenodoMock()
    server.start()
    yield server
    server.stop()


def test_zenodo_snapshot(server: ZenodoMock) -> None:
    run_snakemake(
        "tests/projects/plugins/staging/zenodo-snapshot",
        "staging__upload",
        "--config",
        f"zenodo_mock_url={server.url}/api",
        env=dict(os.environ, ZENODO_TOKEN="test"),
    )


def test_zenodo_restore() -> None:
    run_snakemake(
        "tests/projects/plugins/staging/zenodo-restore",
        "output/a.txt",
        "--config",
        "restore=True",
    )
