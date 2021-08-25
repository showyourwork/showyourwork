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
