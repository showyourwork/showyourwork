from pathlib import Path
import platform


# Get user config, if present
try:
    config
except NameError:
    config = {}

# Verbosity
verbose = str(config.get("verbose", "false")).lower() == "true"

# Recognized figure extensions
figexts = config.get(
    "figexts", ["pdf", "png", "eps", "jpg", "jpeg", "gif", "svg", "tiff"]
)

# Paths to be excluded from the tarball
arxiv_tarball_exclude = config.get(
    "arxiv_tarball_exclude",
    ",".join([
        "**/*.py",
        "**/matplotlibrc",
        "**/.gitignore",
        "**/__pycache__"
    ]),
)

# Install latest version of tectonic from source? May be useful for debugging.
tectonic_latest = str(config.get("tectonic_latest", "false")).lower() == "true"
if platform.system() == "Linux":
    tectonic_os_default = "x86_64-unknown-linux-gnu"
elif platform.system() == "Darwin":
    tectonic_os_default = "x86_64-apple-darwin"
else:
    tectonic_os_default = "x86_64-pc-windows-msvc"
tectonic_os = config.get("tectonic_os", tectonic_os_default)

# Figure dependencies
figure_dependencies = config.get("figure_dependencies", {})
for fd in figure_dependencies:
    full_path = (Path("src") / "figures" / fd).absolute()
    if not full_path.exists():
        raise ValueError("Figure script specified in config file does not exist: {}".format(full_path))