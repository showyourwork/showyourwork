from . import paths, git
from pathlib import Path
import os
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

    # Git info for the repo
    config["git_sha"] = git.get_repo_sha()
    config["git_url"] = git.get_repo_url()
    config["git_branch"] = git.get_repo_branch()
    config["github_actions"] = os.getenv("CI", "false") == "true"
    config["github_runid"] = os.getenv("GITHUB_RUN_ID", "")

    # showyourwork version
    with open(paths.showyourwork / "VERSION", "r") as f:
        config["showyourwork_version"] = f.read().replace("\n", "")

    # Path to the user repo
    config["user_abspath"] = paths.user.as_posix()

    # Path to the workflow
    config["workflow_abspath"] = paths.workflow.as_posix()

    # TeX class name
    config["class_name"] = get_class_name(config["ms_name"])

    # TeX auxiliary files
    config["tex_files_in"] = [
        file.relative_to(paths.user).as_posix()
        for file in (paths.resources / "tex").glob("*")
    ]
    config["tex_files_in"] += [
        file.relative_to(paths.user).as_posix()
        for file in (paths.resources / "classes" / config["class_name"]).glob("*")
    ]
    config["tex_files_out"] = [
        (paths.tex.relative_to(paths.user) / Path(file).name).as_posix()
        for file in config["tex_files_in"]
    ]

    # The main tex file and the compiled pdf
    config["ms_tex"] = (
        paths.tex.relative_to(paths.user) / (config["ms_name"] + ".tex")
    ).as_posix()
    config["ms_pdf"] = config["ms_name"] + ".pdf"

    #
    config["config_json"] = (
        (paths.temp / "config.json").relative_to(paths.user).as_posix()
    )

    config["stylesheet"] = (
        (paths.tex / ".showyourwork.tex").relative_to(paths.user).as_posix()
    )

    config["stylesheet_meta_file"] = (
        (paths.tex / ".showyourwork-metadata.tex").relative_to(paths.user).as_posix()
    )

    # Script extensions (TODO)
    config["script_extensions"] = ["py"]

    # Commands (TODO)
    config["scripts"] = {"py": "python {script}"}

    # Overridden in the `preprocess` rule
    config["tree"] = {"figures": {}}
    config["pdf_dependencies"] = []
    config["labels"] = {}