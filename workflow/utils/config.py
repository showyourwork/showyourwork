from . import paths, git
from pathlib import Path
import snakemake
import re

__all__ = ["parse_config"]


def get_class_name(ms_name):
    """
    Infer the document class for the main TeX file.

    """
    with open(paths.tex / f"{ms_name}.tex", "r") as f:
        lines = f.readlines()
        for line in lines:
            match = re.match("[ \t]*\\\documentclass\[?.*?\]?\{(.*?)\}", line)
            if hasattr(match, "groups"):
                name = match.groups()[0]
                break
        else:
            raise ValueError(f"Unable to determine document class in `{ms_name}.tex`.")
        return name


def parse_config():
    """
    Parse the current config and fill in defaults.

    """
    # Get current config
    config = snakemake.workflow.config

    # -- User settings --

    #: Verbosity
    config["verbose"] = str(config.get("verbose", "false")).lower() == "true"

    #: Manuscript name
    config["ms_name"] = config.get("ms_name", "ms")

    # -- Internal settings --

    # Script extensions (TODO)
    config["script_extensions"] = ["py"]

    # Git info for the repo
    config["git_sha"] = git.get_repo_sha()
    config["git_url"] = git.get_repo_url()
    config["git_branch"] = git.get_repo_branch()

    # Path to the workflow
    config["workflow_path"] = paths.workflow.as_posix()

    # TeX class name
    config["class_name"] = get_class_name(config["ms_name"])

    # TeX auxiliary files
    config["tex_files_in"] = [
        file.as_posix() for file in (paths.resources / "tex").glob("*")
    ]
    config["tex_files_in"] += [
        file.as_posix()
        for file in (paths.resources / "classes" / config["class_name"]).glob("*")
    ]
    config["tex_files_out"] = [
        f"src/tex/{Path(file).name}"
        for file in config["tex_files_in"]
        # if not Path(f"src/tex/{Path(file).name}").exists()
    ]
