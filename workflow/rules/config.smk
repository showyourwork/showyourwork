from pathlib import Path

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
    ]),
)

# Figure dependencies
figure_dependencies = config.get("figure_dependencies", {})
for fd in figure_dependencies:
    full_path = (Path("src") / "figures" / fd).absolute()
    if not full_path.exists():
        raise ValueError("Figure script specified in config file does not exist: {}".format(full_path))