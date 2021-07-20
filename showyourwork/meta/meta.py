from ..showyourwork_version import __version__
from ..constants import *
from ..utils import save_json
from .repo import get_repo_metadata
from .scripts import get_script_metadata, get_script_status
import json


def get_metadata(clobber=True):

    if clobber or not (TEMP / PROJECT / "meta.json").exists():

        # Load the metadata
        repo = get_repo_metadata(clobber=False)
        scripts = get_script_metadata(clobber=False)
        meta = dict(repo=repo)

        # Miscellaneous
        meta["version"] = __version__
        meta["gen_tree"] = False
        meta["graphicspath"] = str(USER / "figures" / "@")[:-1]

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
        save_json(meta, TEMP / PROJECT / "meta.json")

        return meta

    else:

        with open(TEMP / PROJECT / "meta.json", "r") as f:
            meta = json.load(f)

        return meta
