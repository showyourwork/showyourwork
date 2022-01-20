from . import paths, git, exceptions
from pathlib import Path
from collections import OrderedDict, ChainMap
import inspect
import os
import snakemake
import re
import sys

__all__ = ["parse_config", "finalize_config", "get_snakemake_variable", "is_make_clean"]


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


def as_dict(x, depth=0, maxdepth=30):
    """
    Replaces nested instances of OrderedDicts with regular dicts in a dictionary.

    This is useful when parsing a config generated from a YAML file with
    inconsistent use of `-`s.

    """
    if depth == 0 and not x:
        return {}
    elif depth > maxdepth:
        with exceptions.no_traceback():
            # TODO: Our custom excepthook isn't working
            raise exceptions.ConfigError()

    if type(x) is list:
        y = dict(ChainMap(*[dict(xi) for xi in x if type(xi) is OrderedDict]))
        z = [xi for xi in x if type(xi) is not OrderedDict]
        if len(z):
            if y:
                x = [y]
            else:
                x = []
            x.extend(z)
        else:
            x = y
    elif type(x) is OrderedDict:
        x = dict(x)

    if type(x) is dict:
        for key, value in x.items():
            x[key] = as_dict(value, depth + 1)

    return x


def parse_config():
    """
    Parse the current config and fill in defaults.

    This step is only run in the preprocessing stage.

    """
    # Get current config
    config = snakemake.workflow.config

    #
    # -- User settings --
    #

    #: Verbosity
    config["verbose"] = str(config.get("verbose", "false")).lower() == "true"

    #: Debug mode
    config["debug"] = str(config.get("debug", "false")).lower() == "true"

    #: Manuscript name
    config["ms_name"] = config.get("ms_name", "ms")

    #: Custom script execution rules
    config["scripts"] = as_dict(config.get("scripts", {}))
    config["scripts"]["py"] = config["scripts"].get("py", "python {script}")

    #: Custom script dependencies
    config["dependencies"] = as_dict(config.get("dependencies", {}))

    #: Zenodo datasets
    config["zenodo"] = as_dict(config.get("zenodo", {}))
    config["zenodo_sandbox"] = as_dict(config.get("zenodo_sandbox", {}))

    #
    # -- Internal settings --
    #

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

    # The parsed config file
    config["config_json"] = (
        (paths.temp / "config.json").relative_to(paths.user).as_posix()
    )

    # Paths to the TeX stylesheets
    config["stylesheet"] = (
        (paths.tex / ".showyourwork.tex").relative_to(paths.user).as_posix()
    )
    config["stylesheet_meta_file"] = (
        (paths.tex / ".showyourwork-metadata.tex").relative_to(paths.user).as_posix()
    )

    # Script extensions
    config["script_extensions"] = list(config["scripts"].keys())

    # Overridden in the `preprocess` rule
    config["tree"] = {"figures": {}}
    config["labels"] = {}

    # Record additional settings
    finalize_config()


def get_snakemake_variable(name, default=None):
    """
    Infer the value of a variable within snakemake.

    This is extremely hacky.
    """
    for level in inspect.stack():
        value = level.frame.f_locals.get(name, None)
        if value is not None:
            return value
    return default


def is_make_clean():
    """
    Returns True if the current workflow was triggered by `make clean`

    """
    if "--delete-all-output" in sys.argv:
        return True
    else:
        return False


def finalize_config():
    """
    Add some extra settings to the config dict.

    This step is run in both the preprocessing stage and before the main build.
    If we ran it only during preprocessing, passing different command line
    flags to `snakemake` on the next build might have no effect if the
    preprocess workflow is cached & not triggered. The same would be true if
    the user git-committed or changed branches in between builds.

    """
    # Get current config
    config = snakemake.workflow.config

    # Git info for the repo
    config["git_sha"] = git.get_repo_sha()
    config["git_url"] = git.get_repo_url()
    config["git_branch"] = git.get_repo_branch()
    config["github_actions"] = os.getenv("CI", "false") == "true"
    config["github_runid"] = os.getenv("GITHUB_RUN_ID", "")

    # Get the values of some internal snakemake variables
    # so we can mimic snakemake functionality
    config["latency_wait"] = get_snakemake_variable("latency_wait", default=5)
    config["assume_shared_fs"] = get_snakemake_variable(
        "assume_shared_fs", default=True
    )
