import filecmp
import re
import shutil
import subprocess

import jinja2
import yaml
from packaging.version import Version

from .. import exceptions, paths
from ..config import parse_syw_spec
from ..logging import get_logger
from ..subproc import get_stdout

try:
    from yaml import CDumper as Dumper
    from yaml import CLoader as Loader
except ImportError:
    # If LibYAML not installed
    from yaml import Dumper, Loader


# Require this version of conda or greater
MIN_CONDA_VERSION = Version("4.12.0")


def run_in_env(command, **kwargs):
    """Run a command in the isolated showyourwork conda environment.

    This function creates and activates an isolated conda environment
    and executes ``command`` in it, with optional ``kwargs`` passed to
    ``subprocess.run``. The conda environment is specified in
    ``showyourwork/workflow/envs/environment.yml`` and consists primarily
    of the dependencies needed to execute ``Snakemake``, including ``mamba``.
    To this environment we add the specific version of ``showyourwork``
    requested in the article repository's ``showyourwork.yml`` file.

    Note that this ensures that the article is built with the version of
    ``showyourwork`` specified by the workflow, and **not** the version
    the user currently has installed.

    To speed things up on future runs, we cache the environment in the
    temporary folder ``~/.showyourwork/env``.

    Args:
        ``command`` (str): The command to run.
        ``kwargs``: Keyword arguments to pass to subprocess.run.

    Returns:
        subprocess.CompletedProcess: The result of the command.

    Raises:
        ``exceptions.CondaNotFoundError``:
            If conda is not found.
        ``exceptions.CondaVersionError``:
            If an incompatible version of conda installed.
        ``exceptions.ShowyourworkException``:
            If the cwd is not the top level of a showyourwork project.
        ``exceptions.ShowyourworkNotFoundError``:
            If the requested version of showyourwork in the config file cannot be found.
    """
    # Logging
    logger = get_logger()

    # Command to set up conda
    try:
        conda_prefix = get_stdout("conda info --base", shell=True).replace(
            "\n", ""
        )
    except:
        raise exceptions.CondaNotFoundError()
    conda_setup = f". {conda_prefix}/etc/profile.d/conda.sh"

    # Get conda version
    conda_version = get_stdout("conda -V", shell=True)
    try:
        conda_version = Version(
            re.match("^conda (.*?)$", conda_version).groups()[0]
        )
    except:
        raise exceptions.CondaVersionError(MIN_CONDA_VERSION)
    if conda_version < MIN_CONDA_VERSION:
        raise exceptions.CondaVersionError(MIN_CONDA_VERSION, conda_version)

    # Various conda environment files
    syw_envfile = paths.showyourwork().envs / "environment.yml"
    workflow_envfile = paths.user().temp / "environment.yml"
    cached_envfile = paths.user().home_temp / "environment.yml"

    # Infer the `showyourwork` version from the user's config file
    if not (paths.user().repo / "showyourwork.yml").exists():
        raise exceptions.ShowyourworkException(
            "No `showyourwork.yml` config file in current working directory. "
            "Are you running `showyourwork` from within your article's "
            "repository?"
        )
    user_config = yaml.load(
        jinja2.Environment(loader=jinja2.FileSystemLoader(paths.user().repo))
        .get_template("showyourwork.yml")
        .render(),
        Loader=Loader,
    )
    syw_spec = parse_syw_spec(user_config.get("version", None))

    # Copy the showyourwork environment file to a temp location,
    # and add the user's requested showyourwork version as a dependency
    # so we can import it within Snakemake
    with open(syw_envfile, "r") as f:
        syw_env = yaml.load(f, Loader=Loader)
    for dep in syw_env["dependencies"]:
        if type(dep) is dict and "pip" in dep:
            dep["pip"].append(syw_spec)
    with open(workflow_envfile, "w") as f:
        print(yaml.dump(syw_env, Dumper=Dumper), file=f)

    # Set up or update our isolated conda env
    if not paths.user().env.exists():
        # Set up a new env and cache the envfile
        logger.info(
            "Creating a new conda environment in ~/.showyourwork/env..."
        )
        get_stdout(
            f"conda env create --strict-channel-priority -p {paths.user().env} -f {workflow_envfile} -q",
            shell=True,
        )
        shutil.copy(workflow_envfile, cached_envfile)
    else:
        # We'll update the env based on our spec file if the current
        # environment differs (based on checking the cached spec file)
        if cached_envfile.exists():
            cache_hit = filecmp.cmp(
                cached_envfile, workflow_envfile, shallow=False
            )
        else:
            cache_hit = False

        if not cache_hit:
            logger.info("Updating conda environment in ~/.showyourwork/env...")
            get_stdout(
                f"conda env update --strict-channel-priority -p {paths.user().env} -f {workflow_envfile} --prune -q",
                shell=True,
            )
            shutil.copy(workflow_envfile, cached_envfile)

    # Command to activate our environment
    conda_activate = f"{conda_setup} && conda activate {paths.user().env}"

    # Command to get the path to the showyourwork installation.
    # This is used to resolve the path to the internal Snakefile.
    get_syw_path = """SYW_PATH=$(python -c "import showyourwork; from pathlib import Path; print(Path(showyourwork.__file__).parent)")"""

    # Run
    return subprocess.run(
        f"{conda_activate} && {get_syw_path} && {command}",
        shell=True,
        **kwargs,
    )
