from helpers import TemporaryShowyourworkRepository, edit_yaml
from showyourwork.subproc import get_stdout
import pytest


class TestDefault(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo."""

    pass


class TestZenodo(TemporaryShowyourworkRepository):
    """Test a repo that downloads data from Zenodo."""

    pass


class TestZenodoCache(TemporaryShowyourworkRepository):
    """Test the Zenodo caching feature."""

    zenodo_cache = True

    def customize(self):
        """Add the necessary dependencies to the config file."""
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/figures/plot_dataset.py": "src/data/dataset.npz"
            }


class TestZenodoDirCache(TemporaryShowyourworkRepository):
    """Test the Zenodo caching feature for entire directories."""

    zenodo_cache = True

    def customize(self):
        """Add the necessary dependencies to the config file."""
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/figures/plot_dataset.py": "src/data/dataset"
            }