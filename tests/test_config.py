from contextlib import contextmanager
from tempfile import NamedTemporaryFile
from typing import Generator

import pytest
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
