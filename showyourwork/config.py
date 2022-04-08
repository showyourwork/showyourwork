from . import paths, git, exceptions
from .meta import is_make_clean, is_make_main
from pathlib import Path
from collections import OrderedDict, ChainMap
import os
import re

try:
    import snakemake
except ModuleNotFoundError:
    snakemake = None

__all__ = ["parse_config"]


def get_class_name(ms_name):
    """Infer the document class used in the main TeX file.

    Args:
        ms_name (str): The manuscript file name (no path or extension).

    Raises:
        exceptions.UnableToInferClassName: If the name cannot be inferred
            by parsing the TeX file.

    Returns:
        str: The LaTeX class used in the manuscript.
    """
    with open(paths.user().tex / f"{ms_name}.tex", "r") as f:
        lines = f.readlines()
        for line in lines:
            match = re.match("[ \t]*\\\documentclass\[?.*?\]?\{(.*?)\}", line)
            if hasattr(match, "groups"):
                name = match.groups()[0]
                break
        else:
            raise exceptions.UnableToInferClassName(ms_name)
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


    """
    # Get current config
    config = snakemake.workflow.config

    # During the preprocessing stage, we import user settings,
    # set defaults, and record additional internal settings.
    # These get recorded in a JSON file which is loaded as
    # the main config file in the build stage, so all these
    # settings are available in both stages.
    if (
        Path(snakemake.workflow.workflow.main_snakefile).name
        == "preprocess.smk"
    ):

        #
        # -- User settings --
        #

        #: Showyourwork meta-settings
        config["showyourwork"] = config.get("showyourwork", {})
        config["showyourwork"]["version"] = config["showyourwork"].get(
            "version", None
        )
        config["showyourwork"]["cache"] = config["showyourwork"].get(
            "cache", {}
        )
        config["showyourwork"]["cache"]["zenodo"] = config["showyourwork"][
            "cache"
        ].get("zenodo", None)

        #: Verbosity
        config["verbose"] = config.get("verbose", False)

        #: Manuscript name
        config["ms_name"] = config.get("ms_name", "ms")

        #: Custom script execution rules
        config["scripts"] = as_dict(config.get("scripts", {}))
        config["scripts"]["py"] = config["scripts"].get(
            "py", "python {script}"
        )

        #: Custom script dependencies
        config["dependencies"] = as_dict(config.get("dependencies", {}))

        #: Zenodo datasets
        config["zenodo"] = as_dict(config.get("zenodo", {}))
        config["zenodo_sandbox"] = as_dict(config.get("zenodo_sandbox", {}))

        #: Overleaf
        config["overleaf"] = as_dict(config.get("overleaf", {}))

        #
        # -- Internal settings --
        #

        # Path to the user repo
        config["user_abspath"] = paths.user().repo.as_posix()

        # Path to the workflow
        config["workflow_abspath"] = paths.showyourwork().workflow.as_posix()

        # TeX class name
        config["class_name"] = get_class_name(config["ms_name"])

        # TeX auxiliary files
        config["tex_files_in"] = [
            file.as_posix()
            for file in (paths.showyourwork().resources / "tex").glob("*")
        ]
        config["tex_files_in"] += [
            file.as_posix()
            for file in (
                paths.showyourwork().resources
                / "classes"
                / config["class_name"]
            ).glob("*")
        ]
        config["tex_files_out"] = [
            (
                paths.user().tex.relative_to(paths.user().repo)
                / Path(file).name
            ).as_posix()
            for file in config["tex_files_in"]
        ]

        # The main tex file and the compiled pdf
        config["ms_tex"] = (
            paths.user().tex.relative_to(paths.user().repo)
            / (config["ms_name"] + ".tex")
        ).as_posix()
        config["ms_pdf"] = config["ms_name"] + ".pdf"

        # The parsed config file
        config["config_json"] = (
            (paths.user().temp / "config.json")
            .relative_to(paths.user().repo)
            .as_posix()
        )

        # Paths to the TeX stylesheets
        config["stylesheet"] = (
            (paths.user().tex / "showyourwork.tex")
            .relative_to(paths.user().repo)
            .as_posix()
        )
        config["stylesheet_meta_file"] = (
            (paths.user().tex / "showyourwork-metadata.tex")
            .relative_to(paths.user().repo)
            .as_posix()
        )

        # Script extensions
        config["script_extensions"] = list(config["scripts"].keys())

        # Overridden in the `preprocess` rule
        config["tree"] = {"figures": {}}
        config["labels"] = {}
        config["cached_dependencies"] = []

    # The following is run in both the preprocessing stage and the main build.
    # If we ran it only during preprocessing, passing different command line
    # flags to `snakemake` on the next build might have no effect if the
    # preprocess workflow is cached & not triggered. The same would be true if
    # the user git-committed or changed branches in between builds.

    # Git info for the repo
    config["git_sha"] = git.get_repo_sha()
    config["git_url"] = git.get_repo_url()
    config["git_branch"] = git.get_repo_branch()
    config["github_actions"] = os.getenv("CI", "false") == "true"
    config["github_runid"] = os.getenv("GITHUB_RUN_ID", "")