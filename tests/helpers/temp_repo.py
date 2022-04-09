from showyourwork import gitapi
from showyourwork.git import get_repo_sha
from showyourwork.subproc import get_stdout
from showyourwork.zenodo import delete_deposit
from pathlib import Path
import shutil
import asyncio
import pytest
import re
import yaml
import jinja2


# Folder where we'll create our temporary repos
SANDBOX = Path(__file__).absolute().parents[1] / "sandbox"


class TemporaryShowyourworkRepository:
    """
    Base class for showyourwork tests.

    """

    def customize(self):
        """Subclass me to customize the repo."""
        pass

    def create_local(self, zenodo_cache=False, overleaf_id=None):
        """Create the repo locally."""
        # Delete any local repos
        self.delete_local()

        # Parse options
        command = "showyourwork setup"
        options = f"--yes --no-git --showyourwork-version={get_repo_sha()}"
        if not zenodo_cache:
            # Disable zenodo caching
            command = f"ZENODO_TOKEN='' {command}"
        if overleaf_id:
            # Enable overleaf syncing
            options += f"--overleaf={overleaf_id}"

        # Create a new one
        print(
            f"[{self.repo}] Creating local repo `tests/sandbox/{self.repo}`..."
        )
        get_stdout(
            f"{command} {options} showyourwork/{self.repo}",
            cwd=SANDBOX,
            shell=True,
        )

        # Get the Zenodo concept id (if any)
        user_config = yaml.load(
            jinja2.Environment(
                loader=jinja2.FileSystemLoader(SANDBOX / self.repo)
            )
            .get_template("showyourwork.yml")
            .render(),
            Loader=yaml.CLoader,
        )
        self.concept_id = (
            user_config.get("showyourwork", {}).get("cache", {}).get("zenodo")
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

    def setup_git(self):
        """Init the git repo and add + commit all files."""
        print(f"[{self.repo}] Setting up local git repo...")
        get_stdout("git init -q", shell=True, cwd=SANDBOX / self.repo)
        get_stdout("git add .", shell=True, cwd=SANDBOX / self.repo)
        get_stdout(
            "git -c user.name='gh-actions' -c user.email='gh-actions' "
            "commit -q -m 'first commit'",
            shell=True,
            cwd=SANDBOX / self.repo,
        )
        get_stdout("git branch -M main", shell=True, cwd=SANDBOX / self.repo)

    def build_local(self):
        """Run showyourwork locally to build the article."""
        print(f"[{self.repo}] Building the article locally...")
        get_stdout("showyourwork build", shell=True, cwd=SANDBOX / self.repo)

    @pytest.mark.asyncio_cooperative
    async def run_github_action(self, initwait=240, maxtries=10, interval=60):
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
            cwd=SANDBOX / self.repo,
            secrets=[gitapi.get_access_token()],
        )
        print(
            f"[{self.repo}] Waiting {initwait} seconds for workflow "
            f"to finish (1/{maxtries})..."
        )
        await asyncio.sleep(initwait)
        status = "unknown"
        for n in range(maxtries):
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
            elif n < maxtries - 1:
                print(
                    f"[{self.repo}] Waiting {interval} seconds for "
                    f"workflow to finish ({n+2}/{maxtries})..."
                )
                await asyncio.sleep(interval)
        else:
            raise Exception(
                "[{self.repo}] GitHub Actions workflow timed out.\n"
                f"For details, see {url}."
            )

    def delete_zenodo(self):
        """Delete the Zenodo deposit associated with the temp repo."""
        if self.concept_id:
            print(
                f"[{self.repo}] Deleting Zenodo deposit "
                f"with concept id {self.concept_id}..."
            )
            delete_deposit(self.concept_id)

    def delete_remote(self):
        """Delete the remote repo."""
        print(
            f"[{self.repo}] Deleting GitHub repo "
            f"`showyourwork/{self.repo}`..."
        )
        gitapi.delete_repo(self.repo, org="showyourwork")

    def delete_local(self):
        """Delete the local repo."""
        if (SANDBOX / self.repo).exists():
            print(
                f"[{self.repo}] Deleting local repo `tests/sandbox/{self.repo}`..."
            )
            shutil.rmtree(SANDBOX / self.repo)

    @pytest.mark.asyncio_cooperative
    async def test_repo(self):
        """
        Test functionality by creating a new repo, customizing it,
        pushing it to GitHub, and awaiting the article build action to
        complete.

        """
        # Name the repo after the subclass (in hyphenated-snake-case)
        self.repo = (
            re.sub(r"(?<!^)(?=[A-Z])", "_", self.__class__.__name__)
            .lower()
            .replace("_", "-")
        )

        try:

            # Create the repo on GitHub
            self.create_remote()

            # Set up the repo
            self.create_local()

            # Customize the repo
            self.customize()

            # git init, add, and commit
            self.setup_git()

            # Build the article locally
            self.build_local()

            # Push to GitHub to trigger the Actions workflow
            # and wait for the result
            await self.run_github_action()

        except:

            raise

        else:

            # Delete remote repo (only on success)
            self.delete_remote()

            # Delete local repo (only on success)
            self.delete_local()

        finally:

            # Always delete the Zenodo deposit (if created)
            self.delete_zenodo()