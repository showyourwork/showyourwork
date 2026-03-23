import os

import pytest
from helpers import (
    ShowyourworkRepositoryActions,
    TemporaryShowyourworkRepository,
)

from showyourwork import exceptions
from showyourwork.config import edit_yaml
from showyourwork.subproc import get_stdout
from showyourwork.zenodo import Zenodo

pytestmark = pytest.mark.remote


class TestCache(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """
    Test the Zenodo Sandbox caching feature for a rule with a single output.

    """

    # Enable caching
    cache = True
    require_local_build = True

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
            config["run_cache_rules_on_ci"] = True

        # Add the figure environment to the tex file
        self.add_figure_environment()


class TestDirCache(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """
    Test the Zenodo Sandbox caching feature for a rule that outputs an
    entire directory.

    """

    # Enable caching
    cache = True
    require_local_build = True

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
            config["run_cache_rules_on_ci"] = True

        # Add the figure environment to the tex file
        self.add_figure_environment()


if os.getenv("CI", "false") != "true":
    # These tests generate (semi-)permanent deposits on Zenodo Sandbox
    # so it should not be run too often.
    # NOTE: I tried using `pytest.mark.skipif` but they were getting run
    # on CI anyways!

    class TestCacheFreeze(
        TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
    ):
        """Test the Zenodo Sandbox cache freeze feature."""

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
                config["run_cache_rules_on_ci"] = True

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

            # Edit the pipeline script and re-build with seed=1
            # This should download from the cache with seed=1
            self.add_pipeline_script(seed=1)
            self.git_commit()
            self.build_local(env={"SYW_NO_RUN": "true"})

            # Clean so we delete the local cache
            get_stdout("showyourwork clean", cwd=self.cwd, shell=True)

            # Revert the script to `seed=0` and re-build, telling showyourwork
            # to raise an exception if the rule actually gets executed.
            # The expected behavior is for showyourwork to find the frozen version
            # of the cache and use that to bypass the pipeline script run.
            self.add_pipeline_script(seed=0)
            self.git_commit()
            self.build_local(env={"SYW_NO_RUN": "true"})

        def delete_zenodo(self):
            # Override the parent delete_zenodo()
            # We cannot delete the frozen cache
            pass

    class TestCachePublish(
        TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
    ):
        """Test the Zenodo Sandbox cache publish feature."""

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
                config["run_cache_rules_on_ci"] = True

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
            self.build_local(env={"SYW_NO_RUN": "true"})

    class TestCacheReservePublish(
        TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
    ):
        """Test the Zenodo Sandbox cache publish feature."""

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
                config["run_cache_rules_on_ci"] = True

            # Add the figure environment to the tex file
            self.add_figure_environment()

        def check_build(self):
            # Reserve a DOI without publishing anything
            # This normally creates a record *Zenodo*, but we can hack it
            # to use Sandbox so we don't create an actual DOI every time we test!
            get_stdout(
                "SANDBOX_ONLY=true showyourwork cache reserve", cwd=self.cwd, shell=True
            )
            with edit_yaml(self.cwd / "zenodo.yml") as config:
                zenodo_doi = config["cache"].get("main", {}).get("zenodo", None)
            # Make sure the doi was added to the config and that the deposit exists
            assert zenodo_doi is not None
            zenodo = Zenodo(zenodo_doi)
            assert zenodo.check_if_user_is_owner()

            # Make sure that no files have been uploaded
            draft = zenodo._get_draft()
            assert len(draft["files"]) == 0
            assert "notes" not in draft["metadata"]

            # Test that reserving with existing DOI raises an error
            # and leaves the DOI unchanged
            with pytest.raises(
                exceptions.CalledProcessError, match="There is already a DOI"
            ):
                get_stdout(
                    "SANDBOX_ONLY=true showyourwork cache reserve",
                    cwd=self.cwd,
                    shell=True,
                )
            with edit_yaml(self.cwd / "zenodo.yml") as config:
                zenodo_doi_twice = config["cache"].get("main", {}).get("zenodo", None)
            assert zenodo_doi == zenodo_doi_twice

            # Now publish and ensure the same DOI was used and files were uploaded
            get_stdout(
                "SANDBOX_ONLY=true showyourwork cache publish", cwd=self.cwd, shell=True
            )
            with edit_yaml(self.cwd / "zenodo.yml") as config:
                zenodo_doi_pub = config["cache"].get("main", {}).get("zenodo", None)
            assert zenodo_doi_pub == zenodo_doi
            # Files are uploaded
            assert len(draft["files"]) == 0
            assert "notes" not in draft["metadata"]
