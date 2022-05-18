from . import paths, git, exceptions
from pathlib import Path
from collections import OrderedDict, ChainMap
from contextlib import contextmanager
import os
import re
import sys
import jinja2
import yaml

try:
    import snakemake
except ModuleNotFoundError:
    snakemake = None


@contextmanager
def edit_yaml(file):
    """A context used to edit a YAML file in place."""
    if Path(file).exists():
        with open(file, "r") as f:
            contents = yaml.load(f, Loader=yaml.CLoader)
    else:
        contents = {}
    try:
        yield contents
    finally:
        with open(file, "w") as f:
            print(yaml.dump(contents, Dumper=yaml.CDumper), file=f)


def render_config(cwd="."):
    """
    Render any jinja templates in `showyourwork.yml`, combine with
    `zenodo.yml`, and save the processed config to a temporary YAML file.

    This temporary YAML file is then used as the configfile for the Snakemake
    workflow.

    Returns the config as a dictionary.

    """
    # Render the user's config file to a dict
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(cwd))
    config = env.get_template("showyourwork.yml").render()
    config = yaml.safe_load(config)

    # Merge with the zenodo config file, if present
    file = Path(cwd) / "zenodo.yml"
    if file.exists():
        with open(file, "r") as f:
            config.update(yaml.safe_load(f.read()))

    # Save to a temporary YAML file
    with open(paths.user().temp / "showyourwork.yml", "w") as f:
        yaml.dump(config, f)

    return config


def get_run_type():
    """Return the type of the current Snakemake run.

    Options are:

        - ``clean``
        - ``build``
        - ``tarball``
        - ``preprocess``
        - ``other``

    """
    return os.getenv("SNAKEMAKE_RUN_TYPE", "other")


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
    inconsistent use of hyphens.

    """
    if depth == 0 and not x:
        return {}
    elif depth > maxdepth:
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


def parse_overleaf():
    # Get the config
    config = snakemake.workflow.config

    # Make sure `id` is defined
    config["overleaf"]["id"] = config["overleaf"].get("id", None)

    # Make sure `push` and `pull` are defined and they are lists
    config["overleaf"]["push"] = config["overleaf"].get("push", [])
    if config["overleaf"]["push"] is None:
        config["overleaf"]["push"] = []
    elif type(config["overleaf"]["push"]) is not list:
        raise exceptions.ConfigError(
            "Error parsing the config. "
            "The `overleaf.push` field must be a list."
        )
    config["overleaf"]["pull"] = config["overleaf"].get("pull", [])
    if config["overleaf"]["pull"] is None:
        config["overleaf"]["pull"] = []
    elif type(config["overleaf"]["pull"]) is not list:
        raise exceptions.ConfigError(
            "Error parsing the config. "
            "The `overleaf.pull` field must be a list."
        )

    # Ensure all files in `push` and `pull` are in the `src/tex` directory
    for file in config["overleaf"]["push"] + config["overleaf"]["pull"]:
        if not Path(file).resolve().is_relative_to(paths.user().tex):
            raise exceptions.ConfigError(
                "Error parsing the config. "
                "Files specified in `overleaf.push` and `overleaf.pull` must "
                "be located under the `src/tex` directory."
            )

    # Ensure no overlap between `push` and `pull`.
    # User could in principle provide a directory in one
    # and a file within that directory in the other and that would
    # not trigger this error; we'll just have to let them live
    # dangerously!
    push_files = set(
        [
            str(Path(file).resolve().relative_to(paths.user().tex))
            for file in config["overleaf"]["push"]
        ]
    )
    pull_files = set(
        [
            str(Path(file).resolve().relative_to(paths.user().tex))
            for file in config["overleaf"]["pull"]
        ]
    )
    if len(push_files & pull_files):
        raise exceptions.ConfigError(
            "Error parsing the config. "
            "One more more files are listed in both `overleaf.push` and "
            "`overleaf.pull`, which is not supported."
        )


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
    if Path(snakemake.workflow.workflow.main_snakefile).name == "prep.smk":

        #
        # -- User settings --
        #

        #: Showyourwork version
        config["version"] = config.get("version", None)

        #: Verbosity
        config["verbose"] = config.get("verbose", False)

        #: Manuscript name
        config["ms_name"] = config.get("ms_name", "ms")

        #: Custom script execution rules
        # Note that these are strings with placeholder values like `script`
        # inside curly braces; these get formatted later, in `preprocess.py`
        # Note also that since we execute everything from the root of the
        # repo, matplotlib won't by default use the `matplotlibrc` file in
        # the scripts directory, so we need to pass it as the special env
        # var $MATPLOTLIBRC.
        config["scripts"] = as_dict(config.get("scripts", {}))
        config["scripts"]["py"] = config["scripts"].get(
            "py",
            f"MATPLOTLIBRC={paths.user().scripts} " + "python {script}",
        )

        #: Custom script dependencies
        config["dependencies"] = as_dict(config.get("dependencies", {}))

        #: Zenodo datasets
        config["datasets"] = as_dict(config.get("datasets", {}))

        #: Overleaf
        config["overleaf"] = as_dict(config.get("overleaf", {}))
        parse_overleaf()

        #: Latex style customization
        config["style"] = config.get("style", {})
        config["style"]["show_git_sha_or_tag"] = config["style"].get(
            "show_git_sha_or_tag", False
        )

        #: Require inputs to all rules to be present on disk for build to pass?
        config["require_inputs"] = config.get("require_inputs", True)

        #: Allow cached rules to run on CI if there's a cache miss?
        config["run_cache_rules_on_ci"] = config.get(
            "run_cache_rules_on_ci", False
        )

        #
        # -- Internal settings --
        #

        # Cache settings
        config["cache"] = config.get("cache", {})

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
    config["git_tag"] = git.get_repo_tag()
    if config["style"]["show_git_sha_or_tag"]:
        if config["git_tag"] != "":
            config["sha_tag_header"] = f'Git tag: {config["git_tag"]}'
        else:
            # The git default short hash is the first 7 characters:
            config["sha_tag_header"] = f'Git commit: {config["git_sha"][:7]}'
    else:
        config["sha_tag_header"] = ""
    config["cache"][config["git_branch"]] = config["cache"].get(
        config["git_branch"], {}
    )
    config["cache"][config["git_branch"]]["zenodo"] = config["cache"][
        config["git_branch"]
    ].get("zenodo", None)
    config["cache"][config["git_branch"]]["sandbox"] = config["cache"][
        config["git_branch"]
    ].get("sandbox", None)
