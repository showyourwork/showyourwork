from .. import __version__
from .. import paths
from .. import exceptions
from ..logging import get_logger
import subprocess
import shutil
import yaml
import filecmp


def run(command, **kwargs):
    """Run a command in the isolated showyourwork conda environment."""

    # Logging
    logger = get_logger()

    # Command to set up conda
    try:
        conda_prefix = subprocess.check_output("conda info --base", shell=True)
    except subprocess.CalledProcessError:
        # TODO
        raise exceptions.ShowyourworkException()
    conda_prefix = conda_prefix.decode().replace("\n", "")
    conda_setup = f". {conda_prefix}/etc/profile.d/conda.sh"

    # Various conda environment files
    user_envfile = paths.user().repo / "environment.yml"
    syw_envfile = paths.showyourwork().module / "environment.yml"
    workflow_envfile = paths.user().temp / "environment.yml"
    cached_envfile = paths.user().home_temp / "environment.yml"

    # Infer the pip specs for `showyourwork` from the user's env file
    try:
        with open(user_envfile, "r") as f:
            user_env = yaml.load(f, Loader=yaml.CLoader)
    except:
        # TODO
        raise exceptions.ShowyourworkException()
    for dep in user_env["dependencies"]:
        if type(dep) is dict and "pip" in dep:
            pip_deps = dep["pip"]
            break
    else:
        # TODO
        raise exceptions.ShowyourworkException()
    for dep in pip_deps:
        if "showyourwork" in dep:
            syw_spec = dep
            break
    else:
        # TODO
        raise exceptions.ShowyourworkException()

    # Copy the showyourwork environment file to a temp location,
    # and add the user's requested showyourwork version as a dependency
    # so we can import it within Snakemake
    with open(syw_envfile, "r") as f:
        syw_env = yaml.load(f, Loader=yaml.CLoader)
    for dep in syw_env["dependencies"]:
        if type(dep) is dict and "pip" in dep:
            dep["pip"].append(syw_spec)
    with open(workflow_envfile, "w") as f:
        print(yaml.dump(syw_env, Dumper=yaml.CDumper), file=f)

    # Set up or update our isolated conda env
    if not paths.user().env.exists():
        # Set up a new env and cache the envfile
        logger.info(
            "Creating a new conda environment in ~/.showyourwork/env..."
        )
        subprocess.run(
            f"conda env create -p {paths.user().env} -f {workflow_envfile} -q",
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
            subprocess.run(
                f"conda env update -p {paths.user().env} -f {workflow_envfile} --prune -q",
                shell=True,
            )
            shutil.copy(workflow_envfile, cached_envfile)

    # Command to activate our environment
    conda_activate = f"{conda_setup} && conda activate {paths.user().env}"

    # Run
    return subprocess.run(
        f"{conda_activate} && {command}",
        shell=True,
        **kwargs,
    )