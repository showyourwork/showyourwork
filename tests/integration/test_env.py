import os

from helpers import ShowyourworkRepositoryActions, TemporaryShowyourworkRepository


# TODO: Run locally and make sure they both work
# TODO: Since output is captured, the simplest thing is probably to include a
# niche package in the environment and test that it compiles with it and does not w/o
# TODO: Run on CI
# TODO: Then would also be nice to have a way to compile it without conda and not
# get an error
class TestConda(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """Test that showyourwork creates a conda environment by default"""

    def customize(self):
        # Add the script to generate the figure
        self.add_figure_script()

        # Add the figure environment to the tex file
        self.add_figure_environment()

    def check_build(self):
        conda_dir = self.cwd / ".snakemake/conda"
        assert conda_dir.is_dir()
        assert any(conda_dir.iterdir())


class TestNoConda(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """Test that the --no-conda flag works"""

    def customize(self):
        # Add the script to generate the figure
        self.add_figure_script()

        # Add the figure environment to the tex file
        self.add_figure_environment()

        if os.getenv("CI", "false") == "true":
            get_stdout("python -m pip install numpy matplotlib", self.cwd, shell=True)

    def build_local(self):
        super().build_local(args=["--no-conda"])

    def check_build(self):
        conda_dir = self.cwd / ".snakemake/conda"
        assert not any(conda_dir.iterdir())
