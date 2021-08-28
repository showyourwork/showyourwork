# Get user config, if present
try:
    config
except NameError:
    config = {}

# Verbosity
verbose = config.get("verbose", "false").lower() == "true"

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
