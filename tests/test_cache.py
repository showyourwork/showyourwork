from helpers import (
    TemporaryShowyourworkRepository,
    ShowyourworkRepositoryActions,
)
from showyourwork.config import edit_yaml
from showyourwork.subproc import get_stdout
from showyourwork.git import get_commit_message
import pytest
import os


class TestCache(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """Test the Zenodo caching feature."""

    # Enable caching
    cache = True

    def customize(self):
        """Add all necessary files for the build."""
        # Add the pipeline script
        self.add_pipeline_script()

        # Add the Snakefile rule to generate the dataset
        self.add_pipeline_rule()

        # Add the script to generate the figure
        self.add_figure_script(load_data=True)

        # Make the dataset a dependency of the figure
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/test_figure.py": "src/data/test_data.npz"
            }

        # Add the figure environment to the tex file
        self.add_figure_environment()


class TestDirCache(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """Test the Zenodo caching feature for entire directories."""

    # Enable caching
    cache = True

    def customize(self):
        """Add all necessary files for the build."""
        # Add the pipeline script
        self.add_pipeline_script(batch=True)

        # Add the Snakefile rule to generate the dataset
        self.add_pipeline_rule(batch=True)

        # Add the script to generate the figure
        self.add_figure_script(load_data=True, batch=True)

        # Make the dataset a dependency of the figure
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/test_figure.py": "src/data/test_data"
            }

        # Add the figure environment to the tex file
        self.add_figure_environment()


@pytest.mark.skipif(
    os.getenv("CI", "false") == "true",
    reason="test generates a (semi-)permanent deposit on Zenodo Sandbox, "
    "so it should not be run too often",
)
class TestCacheFreeze(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """Test the Zenodo cache freeze feature."""

    # Enable caching
    cache = True

    # No need to test this on CI
    local_build_only = True

    def customize(self):
        """Add all necessary files for the build."""
        # Add the pipeline script
        self.add_pipeline_script(seed=0)

        # Add the Snakefile rule to generate the dataset
        self.add_pipeline_rule()

        # Add the script to generate the figure
        self.add_figure_script(load_data=True)

        # Make the dataset a dependency of the figure
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/test_figure.py": "src/data/test_data.npz"
            }

        # Add the figure environment to the tex file
        self.add_figure_environment()

    def check_build(self):
        # Freeze the cache for `seed=0`, which publishes the cache on Sandbox
        get_stdout("showyourwork cache freeze", cwd=self.cwd, shell=True)

        # Edit the pipeline script and re-build. This will
        # update the latest draft of the Sandbox cache
        self.add_pipeline_script(seed=1)
        self.git_commit()
        self.build_local()

        # Clean so we delete the local cache
        get_stdout("showyourwork clean", cwd=self.cwd, shell=True)

        # Revert the script to `seed=0` and re-build, telling showyourwork
        # to raise an exception if the rule actually gets executed.
        # The expected behavior is for showyourwork to find the frozen version
        # of the cache and use that to bypass the pipeline script run.
        self.add_pipeline_script(seed=0)
        self.git_commit()
        self.build_local(pre="SYW_NO_RUN=true")


@pytest.mark.skipif(
    os.getenv("CI", "false") == "true",
    reason="test generates a (semi-)permanent deposit on Zenodo Sandbox, "
    "so it should not be run too often",
)
class TestCachePublish(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """Test the Zenodo cache publish feature."""

    # Enable caching
    cache = True

    # No need to test this on CI
    local_build_only = True

    def customize(self):
        """Add all necessary files for the build."""
        # Add the pipeline script
        self.add_pipeline_script(seed=0)

        # Add the Snakefile rule to generate the dataset
        self.add_pipeline_rule()

        # Add the script to generate the figure
        self.add_figure_script(load_data=True)

        # Make the dataset a dependency of the figure
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/test_figure.py": "src/data/test_data.npz"
            }

        # Add the figure environment to the tex file
        self.add_figure_environment()

    def check_build(self):
        # Publish the cache for `seed=0`.
        # This normally publishes the cache on *Zenodo*, but we can hack it
        # to use Sandbox so we don't create an actual DOI every time we test!
        get_stdout(
            "SANDBOX_ONLY=true showyourwork cache publish",
            cwd=self.cwd,
            shell=True,
        )

        # Edit the pipeline script and re-build. This will
        # update the latest draft of the Sandbox cache
        self.add_pipeline_script(seed=1)
        self.git_commit()
        self.build_local()

        # Clean so we delete the local cache
        get_stdout("showyourwork clean", cwd=self.cwd, shell=True)

        # Revert the script to `seed=0` and re-build, telling showyourwork
        # to raise an exception if the rule actually gets executed.
        # The expected behavior is for showyourwork to find the frozen version
        # of the cache and use that to bypass the pipeline script run.
        self.add_pipeline_script(seed=0)
        self.git_commit()
        self.build_local(pre="SYW_NO_RUN=true")