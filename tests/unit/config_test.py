import json
from contextlib import contextmanager
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Generator

import pytest
from showyourwork import paths, test_util
from showyourwork.config import ConfigVersionError, ValidationError, load_config


@contextmanager
def temp_config_file(body: str) -> Generator[str, None, None]:
    with NamedTemporaryFile() as f:
        f.write(body.encode("utf-8"))
        f.flush()
        yield f.name


def test_minimal_valid() -> None:
    with temp_config_file("config_version: 2") as f:
        load_config(f)


def test_missing_version() -> None:
    with temp_config_file("") as f:
        with pytest.raises(ConfigVersionError):
            load_config(f)


def test_invalid_version() -> None:
    with temp_config_file("config_version: 1") as f:
        with pytest.raises(ConfigVersionError):
            load_config(f)


def test_invalid_schema() -> None:
    with temp_config_file(
        """
config_version: 2
conda: 1
"""
    ) as f:
        with pytest.raises(ValidationError):
            load_config(f)


def test_snakemake_config() -> None:
    with test_util.temporary_project() as d:
        test_util.run_snakemake(
            str(paths.package_data("showyourwork", "workflow", "Snakefile")),
            ["syw__dump_config", "--config", f"working_directory={d}"],
            cwd=d,
        )
        with open(Path(d) / "config.json", "r") as f:
            data = json.load(f)
        assert data["working_directory"] == str(d)
