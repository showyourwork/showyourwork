from pathlib import Path
import json
import subprocess
import warnings
import re
import shutil
import os

# Read the YAML file
with open(Path(".github") / "workflows" / "showyourwork.yml", "r") as f:
    contents = f.read()

# Replace `current` with the current version (package)
if "showyourwork-version: current" in contents:

    import showyourwork

    version = showyourwork.__version__

    # Add a version file to the repo for the record
    with open(".showyourwork-version", "w") as f:
        print(version, file=f)

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
        "showyourwork-version: current", f"showyourwork-version: {version}"
    )


else:

    version = re.findall("showyourwork-version: (.*?)\n", contents)[0]

    # Add a version file to the repo for the record
    with open(".showyourwork-version", "w") as f:
        print(version, file=f)

# Replace `latest` with the latest version (action)
if "rodluger/showyourwork-action@latest" in contents:

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


# Apply minimal template if user requested it
if "{{cookiecutter.template}}" == "Minimal":
    shutil.move(Path("tex") / "ms_minimal.tex", Path("tex") / "ms.tex")
    os.remove(Path("figures") / "fractals.py")
else:
    os.remove(Path("tex") / "ms_minimal.tex")
