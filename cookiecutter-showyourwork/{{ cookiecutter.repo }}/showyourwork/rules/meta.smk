import subprocess
import os
import re
import json
from pathlib import Path


def get_script_status(script):
    """
    Return an error code corresponding to the git status of a given script.

    """
    # Check if the file exists
    script = Path(script)
    if not script.exists():
        error = ScriptDoesNotExist
    else:
        # Check if the file is version controlled
        try:
            status = (
                subprocess.check_output(
                    ["git", "ls-files", "--error-unmatch", script],
                    cwd=USER,
                    stderr=subprocess.DEVNULL,
                )
                .decode()
                .replace("\n", "")
            )
        except Exception as e:
            # File is not version controlled
            error = ScriptNotVersionControlled
        else:
            # Check if the file has uncommitted changes
            try:
                status = (
                    subprocess.check_output(
                        ["git", "status", "-s", script],
                        cwd=USER,
                        stderr=subprocess.DEVNULL,
                    )
                    .decode()
                    .replace("\n", "")
                )
                if len(status):
                    raise Exception("Uncommitted changes!")
            except Exception as e:
                # Uncommited changes
                error = ScriptHasUncommittedChanges
            else:
                # File is good!
                error = ScriptUpToDate
    return error


def get_showyourwork_version():
    """
    Return the version of `showyourwork` specified in the GitHub workflow.

    """
    try:
        with open(GITHUB / "workflows" / "showyourwork.yml", "r") as f:
            for line in f.readlines():
                match = re.match("[ \t]*showyourwork-version: ([a-f0-9]*)", line)
                if hasattr(match, "groups"):
                    return match.groups()[0]
            else:
                return "main"
    except Exception as e:
        print(e)
        return "main"


localrules: metadata


rule metadata:
    message:
        "Generating article metadata..."
    input:
        [POSIX(file) for file in (GITHUB / "workflows").glob("showyourwork.yml")],
        POSIX(TEMP / "repo.json"),
        POSIX(TEMP / "scripts.json"),
    output:
        POSIX(TEMP / "meta.json")
    run:
        # Load the metadata
        with open(TEMP / "repo.json", "r") as f:
            repo = json.load(f)
        with open(TEMP / "scripts.json", "r") as f:
            scripts = json.load(f)
        meta = dict(repo=repo)

        # Miscellaneous
        meta["version"] = get_showyourwork_version()
        meta["CI"] = os.getenv("CI", "false") == "true"
        meta["run_id"] = os.getenv("GITHUB_RUN_ID", "")

        # Figure metadata
        meta["status"] = ScriptUpToDate
        meta["labels"] = {}
        for label, value in scripts["figures"].items():
            script = value["script"]
            if script != UNKNOWN_SCRIPT:
                status = get_script_status(script)
                meta["labels"]["{}_script".format(label)] = script
                meta["labels"]["{}_status".format(label)] = str(status)
                meta["status"] = max(meta["status"], status)
        meta["status"] = str(meta["status"])

        # Store as JSON
        with open(TEMP / "meta.json", "w") as f:
            print(json.dumps(meta, indent=4), file=f)
