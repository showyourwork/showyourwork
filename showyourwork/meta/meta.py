from ..showyourwork_version import __version__
from ..constants import *
from ..utils import save_json
from .repo import get_repo_metadata
from .scripts import get_script_metadata, get_script_status
import json
from packaging import version
import os
import subprocess


def get_metadata(clobber=True):

    if clobber or not (TEMP / "meta.json").exists():

        # Load the metadata
        repo = get_repo_metadata(clobber=False)
        scripts = get_script_metadata(clobber=False)
        meta = dict(repo=repo)

        # Get the showyourwork version
        showyourwork_version = __version__
        if (
            version.parse(showyourwork_version).is_devrelease
            or version.parse(showyourwork_version).is_prerelease
        ):
            # Use the SHA instead
            try:
                showyourwork_version = (
                    subprocess.check_output(
                        ["git", "rev-parse", "HEAD"],
                        stderr=subprocess.DEVNULL,
                        cwd=ROOT,
                    )
                    .decode()
                    .replace("\n", "")
                )
            except Exception as e:
                showyourwork_version = "main"

        # Miscellaneous
        meta["version"] = showyourwork_version
        meta["gen_tree"] = False
        meta["graphicspath"] = str(USER / "figures" / "@")[:-1]
        meta["CI"] = os.getenv("CI", "false") == "true"
        meta["run_id"] = os.getenv("GITHUB_RUN_ID", "")

        # Figure metadata
        meta["status"] = ScriptUpToDate
        meta["labels"] = {}
        for label, value in scripts["figures"].items():
            script = value["script"]
            status = get_script_status(USER / script)
            meta["labels"]["{}_script".format(label)] = script
            meta["labels"]["{}_status".format(label)] = str(status)
            meta["status"] = max(meta["status"], status)
        meta["status"] = str(meta["status"])

        # Store as JSON
        save_json(meta, TEMP / "meta.json")

        return meta

    else:

        with open(TEMP / "meta.json", "r") as f:
            meta = json.load(f)

        return meta
