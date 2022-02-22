from . import paths, exceptions, logging
from .subproc import run
import subprocess
import os
import shutil
import snakemake
from pathlib import Path


def get_overleaf_credentials(
    overleaf_email="OVERLEAF_EMAIL", overleaf_password="OVERLEAF_PASSWORD"
):
    """
    Return the user's Overleaf email and password, stored in env vars.

    """
    creds = []
    for key in [overleaf_email, overleaf_password]:
        val = os.getenv(key, None)
        if val is None or not len(val):
            # TODO
            raise exceptions.MissingOverleafCredentials()
        else:
            # TODO: Replace other symbols in user's password? Overleaf recommends
            # not using special symbols in the password, which is pretty silly.
            # I bet we can escape the problematic ones here.
            # See https://www.overleaf.com/learn/how-to/Troubleshooting_git_bridge_problems
            val = val.replace("@", "%40")
            creds.append(val)

    return creds


def setup():

    # Logging
    logger = logging.get_logger()

    # Get the Overleaf project id
    project_id = snakemake.workflow.config["overleaf"].get("id", None)
    if not project_id:
        # No Overleaf integration for this project
        return
    else:
        logger.info("Fetching Overleaf repo...")

    # Set up a fresh temp directory
    if paths.overleaf.exists():
        shutil.rmtree(paths.overleaf)
    paths.overleaf.mkdir()

    # Get the credentials & repo url
    overleaf_email, overleaf_password = get_overleaf_credentials()
    url = f"https://{overleaf_email}:{overleaf_password}@git.overleaf.com/{project_id}"

    # Set up a local version of the repo. We don't clone it to avoid
    # storing the url containing the password in .git/config
    run(["git", "init"], cwd=str(paths.overleaf))

    # Pull from the repo (hide secrets)
    run(
        ["git", "pull", url],
        cwd=str(paths.overleaf),
        secrets=[overleaf_email, overleaf_password],
    )


def push_files(files, project_id):

    # Logging
    logger = logging.get_logger()

    # Repo info
    overleaf_email, overleaf_password = get_overleaf_credentials()
    url = f"https://{overleaf_email}:{overleaf_password}@git.overleaf.com/{project_id}"

    # Process each file
    for file in files:

        # Copy it to the local version of the repo
        file = Path(file).resolve()
        remote_file = paths.overleaf / file.relative_to(paths.tex)
        shutil.copy(file, remote_file)

        # git-add the file
        run(
            ["git", "add", remote_file.relative_to(paths.overleaf)],
            cwd=str(paths.overleaf),
        )

    # Commit callback
    def callback(code, stdout, stderr):
        if stdout:
            logger.debug(stdout)
        file_list = " ".join(files)
        if code != 0:
            if "nothing to commit" in stdout:
                logger.info(f"No changes to commit to Overleaf: {file_list}")
            else:
                with exceptions.no_traceback():
                    raise exceptions.CalledProcessError(stderr)
        else:
            logger.info(f"Pushing changes to Overleaf: {file_list}")

    # Commit!
    result = run(
        [
            "git",
            "-c",
            "user.name='showyourwork'",
            "-c",
            "user.email='showyourwork'",
            "commit",
            "-m",
            "automatic showyourwork update",
        ],
        cwd=str(paths.overleaf),
        callback=callback,
    )

    # Push (again being careful about secrets)
    run(
        ["git", "push", url, "master"],
        cwd=str(paths.overleaf),
        secrets=[overleaf_email, overleaf_password],
    )