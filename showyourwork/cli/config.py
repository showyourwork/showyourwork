import json
from pathlib import Path


def get_config():
    """Get current user config."""
    if (Path.home() / ".showyourwork" / "config.json").exists():
        try:
            with open(Path.home() / ".showyourwork" / "config.json", "r") as f:
                config = json.load(f)
        except json.decoder.JSONDecodeError:
            config = {}
    else:
        config = {}
    return config


def write_config(config):
    """Write user config to disk."""
    if not (Path.home() / ".showyourwork" / "config.json").exists():
        (Path.home() / ".showyourwork").mkdir(exist_ok=True)
    with open(Path.home() / ".showyourwork" / "config.json", "w") as f:
        json.dump(config, f)