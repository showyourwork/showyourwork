from contextlib import contextmanager
from tempfile import NamedTemporaryFile
from typing import Generator

import pytest
from showyourwork.config.config import (
    ConfigVersionError,
    ValidationError,
    load_config,
    normalize_keys,
)


def test_normalize_keys() -> None:
    result = normalize_keys(
        {
            "a-b": 1,
            "c-d": {"e-f": 2, "gh": [{"i-j": 3}, {"k-l": 4}]},
        }
    )
    assert result == {
        "a_b": 1,
        "c_d": {"e_f": 2, "gh": [{"i_j": 3}, {"k_l": 4}]},
    }


@contextmanager
def temp_config_file(body: str) -> Generator[str, None, None]:
    with NamedTemporaryFile() as f:
        f.write(body.encode("utf-8"))
        f.flush()
        yield f.name


@pytest.mark.parametrize("sep", ["-", "_"])
def test_minimal_valid(sep: str) -> None:
    with temp_config_file(f"config{sep}version: 2") as f:
        load_config(f)


def test_missing_version() -> None:
    with temp_config_file("") as f:
        with pytest.raises(ConfigVersionError):
            load_config(f)


def test_invalid_version() -> None:
    with temp_config_file("config-version: 1") as f:
        with pytest.raises(ConfigVersionError):
            load_config(f)


def test_invalid_schema() -> None:
    with temp_config_file(
        """
config-version: 2
conda: 1
"""
    ) as f:
        with pytest.raises(ValidationError):
            load_config(f)


def test_document() -> None:
    with temp_config_file(
        """
config-version: 2
document-dependencies:
    - dep0
documents:
    - path: test1
    - test2
    - path: test3
      dependencies:
        - dep1
        - dep2
"""
    ) as f:
        config = load_config(f)
        assert config["documents"] == {
            "test1": ["dep0"],
            "test2": ["dep0"],
            "test3": ["dep0", "dep1", "dep2"],
        }
