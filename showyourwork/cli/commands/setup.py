from ... import paths, exceptions, overleaf, __version__
from ...subproc import get_stdout
from ...logging import get_logger
from ...zenodo import Zenodo
from cookiecutter.main import cookiecutter
from packaging import version
from lastversion import latest
import time
from pathlib import Path
import shutil
import subprocess


def setup(slug, cache, overleaf_id, ssh, showyourwork_version):
    """Set up a new article repo.

    Args:
        slug (str): Repository slug (user/repo).
        cache (bool): If True, enable caching on Zenodo Sandbox.
        overleaf_id (str or NoneType): Overleaf ID of the article.
        ssh (bool): If True, use SSH to clone the repository. Otherwise, use HTTPS.
        showyourwork_version (str): Version of showyourwork to use.

    """
    # Parse the slug
    user, repo = slug.split("/")
    if Path(repo).exists():
        raise exceptions.ShowyourworkException(
            f"Directory already exists: {repo}."
        )
    name = f"@{user}".replace("_", "")

    # Get current stable version
    if showyourwork_version is None:
        if version.parse(__version__).is_devrelease:
            showyourwork_version = str(
                latest("https://pypi.org/project/showyourwork")
            )
        else:
            showyourwork_version = version.parse(__version__).base_version

    # Create a Zenodo deposit draft for this repo
    if cache:
        deposit_sandbox = Zenodo("sandbox", slug=slug, branch="main")
        cache_sandbox_doi = deposit_sandbox.doi
    else:
        deposit_sandbox = None
        cache_sandbox_doi = ""

    # Set up the repo
    cookiecutter(
        str(paths.showyourwork().cookiecutter),
        no_input=True,
        extra_context={
            "user": user,
            "repo": repo,
            "name": name,
            "showyourwork_version": showyourwork_version,
            "cache_sandbox_doi": cache_sandbox_doi,
            "overleaf_id": overleaf_id,
            "year": time.localtime().tm_year,
        },
    )

    # Set up git
    try:
        get_stdout("git init -q", shell=True, cwd=repo)
        get_stdout("git add .", shell=True, cwd=repo)
        get_stdout(
            "git commit -q -m '[showyourwork] first commit'",
            shell=True,
            cwd=repo,
        )
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
        if cache:
            logger.error("Deleting cache deposit...")
            deposit_sandbox.delete()
        logger.error(f"Deleting directory {repo}...")
        shutil.rmtree(repo)
        raise e

    # Set up repository on overleaf
    if overleaf_id is not None:
        try:
            overleaf.setup_remote(overleaf_id, path=Path(repo).absolute())
        except Exception as e:
            logger = get_logger()
            if cache:
                logger.error("Deleting cache deposit...")
                deposit_sandbox.delete()
            logger.error(f"Deleting directory {repo}...")
            shutil.rmtree(repo)
            raise e
