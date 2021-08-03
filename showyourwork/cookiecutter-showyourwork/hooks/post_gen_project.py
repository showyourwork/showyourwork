from pathlib import Path
import json
import subprocess
import warnings
import showyourwork
from packaging.version import parse as parse_version


# Read the YAML file
with open(Path(".github") / "workflows" / "showyourwork.yml", "r") as f:
    contents = f.read()

# Replace `latest` with the current version (package)
version = parse_version(showyourwork.__version__).base_version
try:
    meta = json.loads(
        subprocess.check_output(
            [
                "curl",
                "-s",
                "https://api.github.com/repos/rodluger/showyourwork/releases",
            ]
        ).decode()
    )
    versions = [m["name"] for m in meta]
    if version not in versions:
        raise Exception(
            f"Version `{version}` of `showyourwork` not found on remote."
        )
except Exception as e:
    # Fallback to latest
    warnings.warn(
        f"Version `{version}` of `showyourwork` not found on remote. Falling back to `latest`."
    )
    version = "latest"
contents = contents.replace(
    "showyourwork-version: latest", f"showyourwork-version: {version}"
)

# Replace `latest` with the latest version (action)
try:
    action_version = json.loads(
        subprocess.check_output(
            [
                "curl",
                "-s",
                "https://api.github.com/repos/rodluger/showyourwork-action/releases/latest",
            ]
        ).decode()
    )["name"]
except Exception as e:
    warnings.warn(
        "Unable to determine latest version of `showyourwork-action`."
    )
    action_version = "latest"
contents = contents.replace(
    "rodluger/showyourwork-action@latest",
    f"rodluger/showyourwork-action@{action_version}",
)

# Update the YAML file
with open(Path(".github") / "workflows" / "showyourwork.yml", "w") as f:
    print(contents, file=f)

# Add a VERSION file for the record
with open(".showyourwork-version", "w") as f:
    print(showyourwork.__version__, file=f)
