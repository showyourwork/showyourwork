from showyourwork import gitapi
from showyourwork.git import get_repo_sha
import subprocess
from pathlib import Path
import shutil
import asyncio
import pytest
import re


SANDBOX = Path(__file__).absolute().parents[0] / "sandbox"


class TemporaryShowyourworkRepository:
    def repo_setup(self, zenodo_cache=False, overleaf_id=None):
        """Create the repo locally."""
        # Delete any local repos
        self.repo_cleanup()

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
        print(f"Creating local repo `tests/sandbox/{self.repo}`...")
        res = subprocess.check_output(
            f"{command} {options} showyourwork/{self.repo}",
            cwd=SANDBOX,
            shell=True,
        )

    def repo_cleanup(self):
        """Delete the local repo."""
        if (SANDBOX / self.repo).exists():
            print(f"Deleting local repo `tests/sandbox/{self.repo}`...")
            shutil.rmtree(SANDBOX / self.repo)

    @pytest.mark.asyncio
    async def test_showyourwork(self, initwait=240, maxtries=10, interval=60):
        """
        TODO: Better error handling, logging

        """
        # Name the repo after the subclass (in snake_case)
        self.repo = re.sub(
            r"(?<!^)(?=[A-Z])", "_", self.__class__.__name__
        ).lower()

        try:

            # Create the repo on GitHub
            print(f"Creating GitHub repo `showyourwork/{self.repo}`...")
            gitapi.create_repo(self.repo, org="showyourwork")

            # Set up the repo
            self.repo_setup()

            # Commit and push to GitHub
            print(f"Pushing to `showyourwork/{self.repo}`...")

            # TODO!
            raise NotImplementedError(
                "Use the API to push to a different repo!"
            )

            # Wait for the workflow to finish
            print(
                f"Waiting {initwait} seconds for workflow "
                f"to finish (1/{maxtries})..."
            )
            await asyncio.sleep(initwait)
            status = "unknown"
            for n in range(maxtries):
                (
                    status,
                    conclusion,
                    url,
                ) = gitapi.get_latest_workflow_run_status(
                    self.repo, org="showyourwork"
                )
                if status == "completed":
                    if conclusion == "success":
                        print("Workflow completed successfully.")
                        break
                    else:
                        raise Exception(
                            "GitHub Actions workflow terminated with "
                            f"status {conclusion}.\nFor details, see {url}."
                        )
                elif n < maxtries - 1:
                    print(
                        f"Waiting {interval} seconds for workflow "
                        f"to finish ({n+2}/{maxtries})..."
                    )
                    await asyncio.sleep(interval)
            else:
                raise Exception(
                    "GitHub Actions workflow timed out.\n"
                    f"For details, see {url}."
                )

        except:

            raise

        else:

            # Delete the repo on GitHub on success
            print(f"Deleting GitHub repo `showyourwork/{self.repo}`...")
            gitapi.delete_repo(self.repo, org="showyourwork")

        finally:

            # Clean up the repo
            self.repo_cleanup()


class TestSetupRepo(TemporaryShowyourworkRepository):
    """Test setting up the default repo."""

    pass