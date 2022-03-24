from .. import __version__
from .. import paths
from ..logging import get_logger
import subprocess
import shutil
import packaging.version
import yaml
import filecmp


def run(command, **kwargs):
    """Run a command in the isolated showyourwork conda environment."""

    # Logging
    logger = get_logger()

    # Command to set up conda
    logger.info("Configuring conda...")
    conda_prefix = (
        subprocess.check_output("conda info --base", shell=True)
        .decode()
        .replace("\n", "")
    )
    conda_setup = f". {conda_prefix}/etc/profile.d/conda.sh"

    # Environment variables for the build
    envvars = []

    # Copy the showyourwork environment file to a temp location.
    # If the current version of showyourwork is a released version,
    # add it to the environment file so we can import it within Snakemake.
    # Otherwise (i.e., in dev mode), we'll hack it by recording the path
    # to the package and manually adding it to `sys.path` within Snakemake
    with open(paths.showyourwork().module / "environment.yml", "r") as f:
        env = yaml.load(f, Loader=yaml.CLoader)
    v = packaging.version.parse(__version__)
    if v.is_devrelease or v.is_prerelease or v.is_postrelease:
        # We'll hack it (dev mode)
        envvars.append(f"SHOWYOURWORK_PATH={paths.showyourwork().module}")
    else:
        # Add the exact version to the `environment.yml` spec file
        for dep in env["dependencies"]:
            if type(dep) is dict and "pip" in dep:
                dep["pip"].append(f"showyourwork=={__version__}")
    envfile = paths.user().temp / "environment.yml"
    with open(envfile, "w") as f:
        print(yaml.dump(env, Dumper=yaml.CDumper), file=f)

    # Set up or update our isolated conda env
    cached_envfile = paths.showyourwork().temp / "environment.yml"
    if not paths.showyourwork().env.exists():
        # Set up a new env and cache the envfile
        logger.info(
            "Creating a new conda environment in ~/.showyourwork/env..."
        )
        subprocess.run(
            f"conda env create -p {paths.showyourwork().env} -f {envfile} -q",
            shell=True,
        )
        shutil.copy(envfile, cached_envfile)
    else:
        # We'll update the env based on our spec file if the current
        # environment differs (based on checking the cached spec file)
        if cached_envfile.exists():
            cache_hit = filecmp.cmp(cached_envfile, envfile, shallow=False)
        else:
            cache_hit = False
        if not cache_hit:
            logger.info("Updating conda environment in ~/.showyourwork/env...")
            subprocess.run(
                f"conda env update -p {paths.showyourwork().env} -f {envfile} --prune -q",
                shell=True,
            )
            shutil.copy(envfile, cached_envfile)

    # Command to activate our environment
    conda_activate = (
        f"{conda_setup} && conda activate {paths.showyourwork().env}"
    )

    # Run
    return subprocess.run(
        f"{conda_activate} && {' '.join(envvars)} {command}",
        shell=True,
        **kwargs,
    )