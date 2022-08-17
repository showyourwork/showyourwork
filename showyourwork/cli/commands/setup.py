import shutil
import time
from pathlib import Path

from cookiecutter.main import cookiecutter
from lastversion import latest
from packaging import version

from ... import __version__, exceptions, overleaf, paths
from ...logging import get_logger
from ...subproc import get_stdout
from ...zenodo import Zenodo


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

        # Check if it's an empty git repository
        def callback(code, stdout, stderr):
            if code != 0:
                # `git log` errors if there are no commits;
                # that's what we want
                return
            else:
                # There is some history; let's not re-write it
                raise exceptions.ShowyourworkException(
                    f"Directory already exists: {repo}."
                )

        # Require nothing but a `.git` folder and no commit history
        contents = set(list(Path(repo).glob("*"))) - set([Path(repo) / ".git"])
        if len(contents) == 0:
            get_stdout(
                "git log 2> /dev/null", shell=True, cwd=repo, callback=callback
            )
        else:
            # Refuse to continue
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
        overwrite_if_exists=True,
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

        # Set up the remote if it doesn't exist
        def callback(code, stdout, stderr):
            if code != 0:
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

        get_stdout(
            "git config --get remote.origin.url",
            shell=True,
            cwd=repo,
            callback=callback,
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
