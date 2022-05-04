import showyourwork
from showyourwork import gitapi
from showyourwork.config import render_config
from showyourwork.logging import get_logger
from showyourwork.git import get_repo_sha
from showyourwork.subproc import get_stdout
from showyourwork.zenodo import Zenodo
import logging
from pathlib import Path
import shutil
import asyncio
import pytest
import re
import os


# Folder where we'll create our temporary repos
SANDBOX = Path(__file__).absolute().parents[1] / "sandbox"
RESOURCES = Path(__file__).absolute().parents[1] / "resources"


class TemporaryShowyourworkRepository:
    """
    Base class for showyourwork tests.

    """

    # Local debug mode
    debug = os.getenv("DEBUG", "False").lower() in ("true", "1", "t")

    # Editable class settings
    zenodo_cache = False
    sandbox_cache = False
    overleaf_id = None
    action_wait = 240
    action_max_tries = 10
    action_interval = 60
    use_local_showyourwork = debug
    local_build_only = debug

    # Internal
    _concept_doi = None

    @property
    def repo(self):
        # Name the repo after the subclass (in hyphenated-snake-case)
        return (
            re.sub(r"(?<!^)(?=[A-Z])", "_", self.__class__.__name__)
            .lower()
            .replace("_", "-")
        )

    @property
    def cwd(self):
        return SANDBOX / self.repo

    def startup(self):
        """Subclass me to run things at startup."""
        pass

    def teardown(self):
        """Subclass me to run things at teardown."""
        pass

    def customize(self):
        """Subclass me to customize the repo."""
        pass

    def check_build(self):
        """Subclass me to run post-build local checks."""
        pass

    def copy_files(self):
        """Copy files over from the template repo, if it exists."""
        if (RESOURCES / self.repo).exists():
            for file in (RESOURCES / self.repo).rglob("*"):
                if not file.is_dir():
                    target = self.cwd / file.relative_to(RESOURCES / self.repo)
                    target.parents[0].mkdir(exist_ok=True, parents=True)
                    shutil.copy(file, target)

    def create_local(self):
        """Create the repo locally."""
        # Delete any local repos
        self.delete_local()

        # Parse options
        command = "showyourwork setup"
        if self.use_local_showyourwork:
            version = str(Path(showyourwork.__file__).parents[1])
        else:
            version = get_repo_sha()
        options = f"--quiet --version={version} "
        if self.sandbox_cache:
            # Enable zenodo sandbox caching
            options += "--cache --sandbox"
        elif self.zenodo_cache:
            # Enable zenodo caching
            options += "--cache"
        if self.overleaf_id:
            # Enable overleaf syncing
            options += f"--overleaf={self.overleaf_id}"

        # Provide git name & email
        if os.getenv("CI", "false") == "true":
            get_stdout(
                "git config --global user.name 'gh-actions'",
                cwd=SANDBOX,
                shell=True,
            )
            get_stdout(
                "git config --global user.email 'gh-actions'",
                cwd=SANDBOX,
                shell=True,
            )

        # Create a new one
        print(
            f"[{self.repo}] Creating local repo `tests/sandbox/{self.repo}`..."
        )
        get_stdout(
            f"{command} {options} showyourwork/{self.repo}",
            cwd=SANDBOX,
            shell=True,
        )

        # Get the Zenodo cache concept doi for the main branch (if any)
        self._concept_doi = (
            render_config(cwd=self.cwd)
            .get("showyourwork", {})
            .get("cache", {})
            .get("main")
        )

    def create_remote(self):
        """Create the repo on GitHub."""
        print("")
        print(
            f"[{self.repo}] Creating GitHub repo "
            f"`showyourwork/{self.repo}`..."
        )
        gitapi.create_repo(
            self.repo,
            org="showyourwork",
            description="Temporary test repository for showyourwork",
            private=False,
        )

    def git_commit(self):
        """Add and commit all files in the local repo."""
        get_stdout("git add .", shell=True, cwd=self.cwd)
        get_stdout(
            "git diff-index --quiet HEAD || "
            "git -c user.name='gh-actions' -c user.email='gh-actions' "
            "commit -q -m 'auto commit from showyourwork tests'",
            shell=True,
            cwd=self.cwd,
        )

    def build_local(self):
        """Run showyourwork locally to build the article."""
        print(f"[{self.repo}] Building the article locally...")

        def callback(code, stdout, stderr):
            if code != 0:
                raise Exception(stdout + "\n" + stderr)

        get_stdout(
            "CI=false showyourwork build",
            shell=True,
            cwd=self.cwd,
            callback=callback,
        )

    @pytest.mark.asyncio_cooperative
    async def run_github_action(self):
        """
        Push to the remote and asynchronously wait for the workflow on
        GitHub Actions to finish.
        On success, returns, otherwise raises an Exception.

        """
        print(f"[{self.repo}] Pushing to `showyourwork/{self.repo}`...")
        get_stdout(
            "git push --force https://x-access-token:"
            f"{gitapi.get_access_token()}"
            f"@github.com/showyourwork/{self.repo} main",
            shell=True,
            cwd=self.cwd,
            secrets=[gitapi.get_access_token()],
        )
        print(
            f"[{self.repo}] Waiting {self.action_wait} seconds for workflow "
            f"to finish (1/{self.action_max_tries})..."
        )
        await asyncio.sleep(self.action_wait)
        status = "unknown"
        for n in range(self.action_max_tries):
            (status, conclusion, url,) = gitapi.get_latest_workflow_run_status(
                self.repo, org="showyourwork"
            )
            if status == "completed":
                if conclusion == "success":
                    print(f"[{self.repo}] Workflow completed successfully.")
                    return
                else:
                    raise Exception(
                        "[{self.repo}] GitHub Actions workflow terminated "
                        f"with status {conclusion}.\n"
                        f"For details, see {url}."
                    )
            elif n < self.action_max_tries - 1:
                print(
                    f"[{self.repo}] Waiting {self.action_interval} seconds for "
                    f"workflow to finish ({n+2}/{self.action_max_tries})..."
                )
                await asyncio.sleep(self.action_interval)
        else:
            raise Exception(
                "[{self.repo}] GitHub Actions workflow timed out.\n"
                f"For details, see {url}."
            )

    def delete_zenodo(self):
        """Delete the Zenodo deposit associated with the temp repo."""
        if self._concept_doi:
            print(
                f"[{self.repo}] Deleting Zenodo deposit "
                f"with concept DOI {self._concept_doi}..."
            )
            Zenodo(self._concept_doi).delete()

    def delete_remote(self):
        """Delete the remote repo."""
        print(
            f"[{self.repo}] Deleting GitHub repo "
            f"`showyourwork/{self.repo}`..."
        )
        gitapi.delete_repo(self.repo, org="showyourwork")

    def delete_local(self):
        """Delete the local repo."""
        if (self.cwd).exists():
            print(
                f"[{self.repo}] Deleting local repo `tests/sandbox/{self.repo}`..."
            )
            shutil.rmtree(self.cwd)

    def disable_logging(self):
        """Disable showyourwork screen output."""
        for handler in get_logger().handlers:
            if isinstance(handler, logging.StreamHandler):
                handler.setLevel(logging.ERROR)

    @pytest.mark.asyncio_cooperative
    async def test_repo(self):
        """
        Test functionality by creating a new repo, customizing it,
        pushing it to GitHub, and awaiting the article build action to
        complete.

        """
        try:

            # Disable screen logging info from showyourwork
            self.disable_logging()

            # Always run this first
            self.startup()

            # Create the repo on GitHub
            if not self.local_build_only:
                self.create_remote()

            # Set up the repo
            self.create_local()

            # Copy files from the template
            self.copy_files()

            # Customize the repo
            self.customize()

            # Commit changes
            self.git_commit()

            # Build the article locally
            self.build_local()

            # Run local checks
            self.check_build()

            # Push to GitHub to trigger the Actions workflow
            # and wait for the result
            if not self.local_build_only:
                await self.run_github_action()

        except:

            raise

        else:

            # Delete remote repo (only on success)
            if not self.local_build_only:
                self.delete_remote()

            # Delete local repo (only on success)
            self.delete_local()

        finally:

            # Always delete the Zenodo deposit (if created)
            self.delete_zenodo()

            # Always run this last
            self.teardown()