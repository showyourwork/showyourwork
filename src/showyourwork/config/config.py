from pathlib import Path
from typing import Any, Dict, Union

import yaml
from importlib.resources import files
from jsonschema import validate, ValidationError as JSONSchemaValidationError

from showyourwork.version import __version__

PathLike = Union[str, Path]


class ConfigVersionError(Exception):
    def __init__(self, config: Dict[str, Any], required_version: int = 2):
        if "config-version" not in config:
            super().__init__(
                "The expected version of the showyourwork configuration must be "
                "specified using the 'config-version' key."
            )
        else:
            super().__init__(
                f"Version {__version__} of showyourwork requires configuration with "
                f"version {required_version} of the specification, but the config file "
                f"specifies the version as {config['config-version']}"
            )


class ValidationError(Exception):
    pass


def parse_config(file: PathLike, required_version: int = 2) -> Dict[str, Any]:
    with open(file, "r") as f:
        config = yaml.safe_load(f)

    if config is None:
        config = {}

    if config.get("config-version", None) != required_version:
        raise ConfigVersionError(config, required_version=required_version)

    with files("showyourwork.config").joinpath("config.schema.yaml").open("r") as f:
        schema = yaml.safe_load(f)

    try:
        validate(config, schema)
    except JSONSchemaValidationError as e:
        raise ValidationError(
            "The configuration file is invalid; schema validation failed with the "
            f"following error:\n\n{e}"
        ) from e

    return config
