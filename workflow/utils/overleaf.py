from . import paths, exceptions, logging
from .subproc import run
import subprocess
import os
import shutil
import snakemake
from pathlib import Path
from urllib.parse import quote


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
            # Replace special characters in the credentials
            val = quote(val, safe="")
            creds.append(val)

    return creds


def clone(project_id):

    # Logging
    logger = logging.get_logger()
    logger.info("Cloning Overleaf repo...")

    # Set up a fresh temp directory
    if paths.overleaf.exists():
        shutil.rmtree(paths.overleaf)
    paths.overleaf.mkdir()

    # Get the credentials & repo url
    overleaf_email, overleaf_password = get_overleaf_credentials()
    url = f"https://{overleaf_email}:{overleaf_password}@git.overleaf.com/{project_id}"

    # Set up a local version of the repo. We don't actually _clone_ it to avoid
    # storing the url containing the password in .git/config
    run(["git", "init"], cwd=str(paths.overleaf))

    # Pull from the repo (hide secrets)
    run(
        ["git", "pull", url],
        cwd=str(paths.overleaf),
        secrets=[overleaf_email, overleaf_password],
    )


def push_files(files, project_id):

    # Setup logging
    logger = logging.get_logger()

    # Clone the repo
    clone(project_id)

    # Process each file
    skip = []
    for file in files:

        # Copy it to the local version of the repo
        if not Path(file).exists():
            skip.append(file)
            continue
        file = Path(file).resolve()
        remote_file = paths.overleaf / file.relative_to(paths.tex)
        if file.is_dir():
            if remote_file.exists():
                shutil.rmtree(remote_file)
            shutil.copytree(file, remote_file)
        else:
            shutil.copy(file, remote_file)

        # git-add the file
        run(
            ["git", "add", remote_file.relative_to(paths.overleaf)],
            cwd=str(paths.overleaf),
        )

    # Remove missing files from the list
    if skip:
        skip_list = " ".join([str(s) for s in skip])
        logger.warn(f"Skipping missing file(s): {skip_list}")
        files = list(set(files) - set(skip))

    # Commit callback
    def callback(code, stdout, stderr):
        if stdout:
            logger.debug(stdout)
        file_list = " ".join([str(s) for s in files])
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
    overleaf_email, overleaf_password = get_overleaf_credentials()
    url = f"https://{overleaf_email}:{overleaf_password}@git.overleaf.com/{project_id}"
    run(
        ["git", "push", url, "master"],
        cwd=str(paths.overleaf),
        secrets=[overleaf_email, overleaf_password],
    )


def pull_files(files, project_id):

    # Clone the repo
    clone(project_id)

    # Copy over the files
    file_list = " ".join(files)
    logger = logging.get_logger()
    logger.info(f"Pulling changes from Overleaf: {file_list}")
    logger.warn("Run `git status` to see what changed. Don't forget to commit!")
    for file in files:
        file = Path(file).absolute()
        remote_file = (paths.overleaf / file.relative_to(paths.tex)).resolve()
        if not remote_file.exists():
            raise exceptions.MissingOverleafFile(
                remote_file.relative_to(paths.overleaf)
            )
        if remote_file.is_dir():
            if file.exists():
                shutil.rmtree(file)
            shutil.copytree(remote_file, file)
        else:
            shutil.copy(remote_file, file)