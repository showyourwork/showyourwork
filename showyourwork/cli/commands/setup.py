from ... import paths, zenodo, exceptions, overleaf, __version__
from ...subproc import get_stdout
from ...logging import get_logger
from cookiecutter.main import cookiecutter
from packaging import version
import time
from pathlib import Path
import shutil
import subprocess


def setup(slug, overleaf_id, ssh, showyourwork_version):
    """Set up a new article repo."""
    # Parse the slug
    user, repo = slug.split("/")
    if Path(repo).exists():
        raise exceptions.ShowyourworkException(
            f"Directory already exists: {repo}."
        )
    name = f"@{user}"

    # Get current stable version
    if showyourwork_version is None:
        if version.parse(__version__).is_devrelease:
            showyourwork_version = "latest"
        else:
            showyourwork_version = version.parse(__version__).base_version

    # Create a Zenodo deposit draft for this repo if the user
    # set the ZENODO_TOKEN environment variable
    if zenodo.get_access_token(error_if_missing=False):
        zenodo_cache_concept_id = zenodo.create_deposit(
            f"Data for {slug} [main]", zenodo_url="zenodo.org"
        )
    else:
        zenodo_cache_concept_id = ""

    # Set up the repo
    cookiecutter(
        str(paths.showyourwork().cookiecutter),
        no_input=True,
        extra_context={
            "user": user,
            "repo": repo,
            "name": name,
            "showyourwork_version": showyourwork_version,
            "zenodo_cache_concept_id": zenodo_cache_concept_id,
            "overleaf_id": overleaf_id,
            "year": time.localtime().tm_year,
        },
    )

    # Set up git
    try:
        get_stdout("git init -q", shell=True, cwd=repo)
        get_stdout("git add .", shell=True, cwd=repo)
        get_stdout("git commit -q -m 'first commit'", shell=True, cwd=repo)
        get_stdout("git branch -M main", shell=True, cwd=repo)
        if ssh:
            get_stdout(
                f"git remote add origin git@github.com:{user}/{repo}.git",
                shell=True,
                cwd=repo,
            )
        else:
            get_stdout(
                f"git remote add origin https://github.com/{user}/{repo}.git",
                shell=True,
                cwd=repo,
            )
    except Exception as e:
        logger = get_logger()
        if zenodo_cache_concept_id:
            logger.error("Deleting Zenodo deposit...")
            zenodo.delete_deposit(zenodo_cache_concept_id)
        logger.error(f"Deleting directory {repo}...")
        shutil.rmtree(repo)
        raise e

    # Set up repository on overleaf
    if overleaf_id is not None:
        try:
            overleaf.setup_remote(overleaf_id, path=Path(repo).absolute())
        except Exception as e:
            logger = get_logger()
            if zenodo_cache_concept_id:
                logger.error("Deleting Zenodo deposit...")
                zenodo.delete_deposit(zenodo_cache_concept_id)
            logger.error(f"Deleting directory {repo}...")
            shutil.rmtree(repo)
            raise e
