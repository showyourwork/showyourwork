from tempfile import NamedTemporaryFile

import pytest

from showyourwork.config import parse_config, ConfigVersionError, ValidationError


def test_minimal_valid() -> None:
    with NamedTemporaryFile() as f:
        f.write(
            """
config-version: 2
        """.encode(
                "utf-8"
            )
        )
        f.flush()
        parse_config(f.name)


def test_missing_version() -> None:
    with NamedTemporaryFile() as f:
        f.write("".encode("utf-8"))
        f.flush()
        with pytest.raises(ConfigVersionError):
            parse_config(f.name)


def test_invalid_version() -> None:
    with NamedTemporaryFile() as f:
        f.write("config-version: 1".encode("utf-8"))
        f.flush()
        with pytest.raises(ConfigVersionError):
            parse_config(f.name)


def test_invalid_schema() -> None:
    with NamedTemporaryFile() as f:
        f.write(
            """
config-version: 2
conda: 1
""".encode(
                "utf-8"
            )
        )
        f.flush()
        with pytest.raises(ValidationError):
            parse_config(f.name)
