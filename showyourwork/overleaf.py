from . import paths, exceptions, logging
from .subproc import run
import subprocess
import os
import shutil
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
            raise exceptions.MissingOverleafCredentials()
        else:
            # Replace special characters in the credentials
            val = quote(val, safe="")
            creds.append(val)

    return creds


def clone(project_id):

    # Logging
    logger = logging.get_logger()
    logger.info("Fetching Overleaf repo...")

    # Set up a fresh temp directory
    if paths.user().overleaf.exists():
        shutil.rmtree(paths.user().overleaf)
    paths.user().overleaf.mkdir()

    # Get the credentials & repo url
    overleaf_email, overleaf_password = get_overleaf_credentials()
    url = f"https://{overleaf_email}:{overleaf_password}@git.overleaf.com/{project_id}"

    # Set up a local version of the repo. We don't actually _clone_ it to avoid
    # storing the url containing the password in .git/config
    run(["git", "init"], cwd=str(paths.user().overleaf))

    # Pull from the repo (hide secrets)
    def callback(code, stdout, stderr):
        if stdout:
            logger.debug(stdout)
        if code != 0:
            if "Authentication failed" in stderr:
                raise exceptions.OverleafAuthenticationError()
            else:
                with exceptions.no_traceback():
                    raise exceptions.CalledProcessError(stdout + "\n" + stderr)

    run(
        ["git", "pull", url],
        cwd=str(paths.user().overleaf),
        secrets=[overleaf_email, overleaf_password],
        callback=callback,
    )


def push_files(files, project_id):

    # Disable if user didn't specify an id or if there are no files
    if not project_id or not files:
        return

    # Setup logging
    logger = logging.get_logger()

    # Clone the repo
    try:
        clone(project_id)
    except (
        exceptions.MissingOverleafCredentials,
        exceptions.OverleafAuthenticationError,
    ):
        # Not fatal!
        return

    # Process each file
    skip = []
    for file in files:

        # Copy it to the local version of the repo
        if not Path(file).exists():
            skip.append(file)
            continue
        file = Path(file).resolve()
        remote_file = paths.user().overleaf / file.relative_to(
            paths.user().tex
        )
        if file.is_dir():
            if remote_file.exists():
                shutil.rmtree(remote_file)
            shutil.copytree(file, remote_file)
        else:
            shutil.copy(file, remote_file)

        # git-add the file
        run(
            ["git", "add", remote_file.relative_to(paths.user().overleaf)],
            cwd=str(paths.user().overleaf),
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
            if (
                "Your branch is up to date" in stdout + stderr
                or "nothing to commit" in stdout + stderr
                or "nothing added to commit" in stdout + stderr
            ):
                logger.warn(f"No changes to commit to Overleaf: {file_list}")
            else:
                raise exceptions.CalledProcessError(stdout + "\n" + stderr)
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
        cwd=str(paths.user().overleaf),
        callback=callback,
    )

    # Push (again being careful about secrets)
    overleaf_email, overleaf_password = get_overleaf_credentials()
    url = f"https://{overleaf_email}:{overleaf_password}@git.overleaf.com/{project_id}"
    run(
        ["git", "push", url, "master"],
        cwd=str(paths.user().overleaf),
        secrets=[overleaf_email, overleaf_password],
    )


def pull_files(files, project_id, auto_commit=False):

    # Disable if user didn't specify an id or if there are no files
    if not project_id or not files:
        return

    # Setup logging
    logger = logging.get_logger()

    # Clone the repo
    try:
        clone(project_id)
    except (
        exceptions.MissingOverleafCredentials,
        exceptions.OverleafAuthenticationError,
    ):
        # Not fatal!
        return

    # Copy over the files
    file_list = " ".join(files)
    logger.info(f"Pulling changes from Overleaf: {file_list}")
    for file in files:
        file = Path(file).absolute()
        remote_file = (
            paths.user().overleaf / file.relative_to(paths.user().tex)
        ).resolve()
        if not remote_file.exists():
            # Non-fatal
            logger.error(
                f"File not found on Overleaf: {remote_file.relative_to(paths.user().overleaf)}"
            )
            continue

        if remote_file.is_dir():
            if file.exists():
                shutil.rmtree(file)
            shutil.copytree(remote_file, file)
        else:
            # Only copy if the files actually differ
            def callback(code, stdout, stderr):
                if code != 0:
                    shutil.copy(remote_file, file)
                    if auto_commit:
                        run(
                            [
                                "git",
                                "add",
                                file.relative_to(paths.user().repo),
                            ],
                            cwd=paths.user().repo,
                        )

            run(["diff", remote_file, file], callback=callback)

    if auto_commit:

        def callback(code, stdout, stderr):
            if code == 0:
                # Push the changes
                def callback(code, stdout, stderr):
                    if code != 0:
                        # Non-fatal
                        logger.error(
                            "Error pushing changes to the remote.\n" + stderr
                        )

                run(["git", "push"], cwd=paths.user().repo, callback=callback)
                logger.info("Changes committed and pushed to remote.")
            else:
                if (
                    "Your branch is up to date" in stdout + stderr
                    or "nothing to commit" in stdout + stderr
                    or "nothing added to commit" in stdout + stderr
                ):
                    logger.warn(f"No changes to commit to the repo.")
                else:
                    raise exceptions.CalledProcessError(stdout + "\n" + stderr)

        run(
            [
                "git",
                "-c",
                "user.name='showyourwork'",
                "-c",
                "user.email='showyourwork'",
                "commit",
                "-m",
                "automatic showyourwork Overleaf update",
            ],
            cwd=paths.user().repo,
            callback=callback,
        )

    else:

        logger.warn(
            "Run `git status` to see what changed. Don't forget to commit!"
        )
