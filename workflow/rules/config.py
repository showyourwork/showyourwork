"""
Process the user config settings.

"""
from pathlib import Path
import platform
import os
from sphinx_mock import *


__all__ = ["config"]


# Get user config, if present
try:
    config
except NameError:
    #: Dictionary containing all the user config settings
    config = {}

#: Verbosity
config["verbose"] = str(config.get("verbose", "false")).lower() == "true"

#: Recognized figure extensions
config["figexts"] = config.get(
    "figexts", ["pdf", "png", "eps", "jpg", "jpeg", "gif", "svg", "tiff"]
)

#: Paths to be excluded from the tarball
config["arxiv_tarball_exclude"] = config.get(
    "arxiv_tarball_exclude",
    [],
)

#: Install latest version of tectonic from source? May be useful for debugging.
config["tectonic_latest"] = (
    str(config.get("tectonic_latest", "false")).lower() == "true"
)
if config["tectonic_latest"]:
    config["tectonic_cmd"] = ".showyourwork/tectonic"
else:
    config["tectonic_cmd"] = "tectonic"
if platform.system() == "Linux":
    tectonic_os_default = "x86_64-unknown-linux-gnu"
elif platform.system() == "Darwin":
    tectonic_os_default = "x86_64-apple-darwin"
else:
    tectonic_os_default = "x86_64-pc-windows-msvc"
config["tectonic_os"] = config.get("tectonic_os", tectonic_os_default)


#: Figure dependencies
config["figure_dependencies"] = config.get("figure_dependencies", {})
for fd in config["figure_dependencies"]:
    full_path = (Path("src") / "figures" / fd).absolute()
    if not full_path.exists():
        raise ValueError(
            "Figure script specified in config file does not exist: {}".format(
                full_path
            )
        )


#: Are we on GitHub Actions?
config["CI"] = os.getenv("CI", "false") == "true"