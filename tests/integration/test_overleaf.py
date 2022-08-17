import time

import pytest
from helpers import (
    ShowyourworkRepositoryActions,
    TemporaryShowyourworkRepository,
)

from showyourwork import exceptions, overleaf
from showyourwork.config import edit_yaml
from showyourwork.logging import get_logger
from showyourwork.subproc import get_stdout

pytestmark = pytest.mark.remote


class TestOverleaf(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """Test a repo that integrates with an Overleaf project."""

    # No need to test this on CI
    local_build_only = True
    overleaf_id = "6262c032aae5421d6d945acf"

    # Overeleaf rate limit error re-try settings
    auth_retries = 1
    auth_sleep = 60

    def startup(self):
        """Wipe the Overleaf remote to start fresh."""
        for n in range(self.auth_retries):
            try:
                overleaf.wipe_remote(self.overleaf_id)
            except exceptions.OverleafRateLimitExceeded:
                get_logger().warn(
                    f"Overleaf authentication failed. Re-trying in {self.auth_sleep} seconds..."
                )
                time.sleep(self.auth_sleep)
            else:
                break

    def customize(self):
        """Customize the build to test all Overleaf functionality."""
        # Add a figure script
        self.add_figure_script()

        # Add Overleaf settings to the config
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["overleaf"]["push"] = ["src/tex/figures"]
            config["overleaf"]["pull"] = ["src/tex/ms.tex"]

        # Mock changes on the Overleaf remote by modifying the local
        # manuscript (adding a figure environment),
        # pushing it to Overleaf, and restoring the original local version
        self.add_figure_environment()
        for n in range(self.auth_retries):
            try:
                overleaf.push_files(
                    [self.cwd / "src" / "tex" / "ms.tex"],
                    self.overleaf_id,
                    path=self.cwd,
                )
            except exceptions.OverleafRateLimitExceeded:
                get_logger().warn(
                    f"Overleaf authentication failed. Re-trying in {self.auth_sleep} seconds..."
                )
                time.sleep(self.auth_sleep)
            else:
                break

        get_stdout("git checkout -- src/tex/ms.tex", shell=True, cwd=self.cwd)

    def check_build(self):
        """Run several post-build checks.

        First, we check if the figure PDF exists. It should only exist if we
        successfully synced the Overleaf manuscript to the local repo, as
        that contains the figure environment defining the script that produces
        `test_figure.pdf`.

        Here we also check that the generated figure was pushed to the
        Overleaf project successfully.

        We then make a change to the local version of the manuscript
        and check that we get an error alerting us to the fact that those
        changes would get overwritten; we then commit the changes and run the
        same check, which should throw a similar error. Finally, we amend the
        commit message to include the magical [showyourwork] stamp, which
        is how we tell showyourwork the file has been manually synced to
        Overleaf (see the docs).
        """
        # Check that the generated figure is present locally
        figure = self.cwd / "src" / "tex" / "figures" / "test_figure.pdf"
        assert figure.exists()

        # Check that the figure is present on the remote
        for n in range(self.auth_retries):
            try:
                overleaf.pull_files(
                    [figure],
                    self.overleaf_id,
                    path=self.cwd,
                    error_if_missing=True,
                )
            except exceptions.OverleafRateLimitExceeded:
                get_logger().warn(
                    f"Overleaf authentication failed. Re-trying in {self.auth_sleep} seconds..."
                )
                time.sleep(self.auth_sleep)
            else:
                break

        # Check that an exception is raised if we try to overwrite a file
        # with uncommitted changes
        ms = self.cwd / "src" / "tex" / "ms.tex"
        with open(ms, "r") as f:
            ms_orig = f.read()
        with open(ms, "w") as f:
            f.write(r"% dummy comment\n" + ms_orig)
        for n in range(self.auth_retries):
            try:
                overleaf.pull_files(
                    [ms],
                    self.overleaf_id,
                    path=self.cwd,
                    error_if_local_changes=True,
                )
            except exceptions.OverleafRateLimitExceeded:
                get_logger().warn(
                    f"Overleaf authentication failed. Re-trying in {self.auth_sleep} seconds..."
                )
                time.sleep(self.auth_sleep)
            except exceptions.OverleafError as e:
                break
            else:
                raise Exception("Failed to raise exception!")

        # Commit the changes and check that the exception is still raised
        get_stdout(
            f"git add -f {ms} && git commit -m 'changing ms.tex locally'",
            cwd=self.cwd,
            shell=True,
        )
        for n in range(self.auth_retries):
            try:
                overleaf.pull_files(
                    [ms],
                    self.overleaf_id,
                    path=self.cwd,
                    error_if_local_changes=True,
                )
            except exceptions.OverleafRateLimitExceeded:
                get_logger().warn(
                    f"Overleaf authentication failed. Re-trying in {self.auth_sleep} seconds..."
                )
                time.sleep(self.auth_sleep)
            except exceptions.OverleafError as e:
                break
            else:
                raise Exception("Failed to raise exception!")

        # Amend the commit message with the magical `[showyourwork]` label
        # and check that the merge works
        get_stdout(
            f"git commit --amend -m '[showyourwork] changing ms.tex locally'",
            cwd=self.cwd,
            shell=True,
        )
        for n in range(self.auth_retries):
            try:
                overleaf.pull_files(
                    [ms],
                    self.overleaf_id,
                    path=self.cwd,
                    error_if_local_changes=True,
                )
            except exceptions.OverleafRateLimitExceeded:
                get_logger().warn(
                    f"Overleaf authentication failed. Re-trying in {self.auth_sleep} seconds..."
                )
                time.sleep(self.auth_sleep)
            else:
                break
