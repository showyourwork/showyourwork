from helpers import (
    TemporaryShowyourworkRepository,
    ShowyourworkRepositoryActions,
)
from showyourwork.config import edit_yaml


class TestSandboxCache(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """Test the Zenodo caching feature."""

    # Enable Sandbox caching
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


class TestSandboxDirCache(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """Test the Zenodo caching feature for entire directories."""

    # Enable Sandbox caching
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