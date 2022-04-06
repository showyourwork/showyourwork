from ... import paths, zenodo, exceptions, __version__
from cookiecutter.main import cookiecutter
from packaging import version
import time
from pathlib import Path
import subprocess


def git_init_commit_and_push(
    user, repo, cwd, ssh=False, authenticate_as_showyourwork=False
):
    subprocess.run("git init -q", shell=True, cwd=cwd)
    if authenticate_as_showyourwork:
        subprocess.run(
            "git config user.name 'showyourwork'", shell=True, cwd=cwd
        )
        subprocess.run(
            "git config user.email 'showyourwork'", shell=True, cwd=cwd
        )
    subprocess.run("git add .", shell=True, cwd=cwd)
    subprocess.run("git commit -q -m 'first commit'", shell=True, cwd=cwd)
    subprocess.run("git branch -M main", shell=True, cwd=cwd)
    if ssh:
        subprocess.run(
            f"git remote add origin git@github.com:{user}/{repo}.git",
            shell=True,
            cwd=cwd,
        )
    else:
        subprocess.run(
            f"git remote add origin https://github.com/{user}/{repo}.git",
            shell=True,
            cwd=cwd,
        )
    res = subprocess.run("git push -q -u origin main", shell=True, cwd=cwd)
    if res.returncode > 0:
        with exceptions.no_traceback():
            raise exceptions.ShowyourworkException(
                f"Unable to push to GitHub. Did you forget to create the remote repo?"
            )


def setup(slug, overleaf, ssh, no_git, showyourwork_version):
    """Set up a new article repo."""
    # Parse the slug
    user, repo = slug.split("/")
    if Path(repo).exists():
        with exceptions.no_traceback():
            raise exceptions.ShowyourworkException(
                f"Directory already exists: {repo}."
            )
    name = f"@{user}"

    # Get current stable version
    if showyourwork_version is None:
        showyourwork_version = version.parse(__version__).base_version

    # Create a Zenodo deposit draft for this repo if the user
    # set the ZENODO_TOKEN environment variable
    if zenodo.get_access_token(error_if_missing=False):
        zenodo_cache_concept_id = zenodo.create_deposit(
            f"Data for {slug} [main]"
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
            "overleaf_id": overleaf,
            "year": time.localtime().tm_year,
        },
    )

    # Set up git
    if not no_git:
        git_init_commit_and_push(user, repo, repo, ssh=ssh)