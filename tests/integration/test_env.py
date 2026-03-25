import os

from helpers import ShowyourworkRepositoryActions, TemporaryShowyourworkRepository

from showyourwork.subproc import get_stdout


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
            get_stdout(
                "/home/runner/micromamba-bin/micromamba install -y tectonic=0.14.1",
                cwd=self.cwd,
                shell=True,
            )
            get_stdout(
                "python -m pip install numpy matplotlib", cwd=self.cwd, shell=True
            )

    def build_local(self):
        super().build_local(args=["--no-conda"])

    def check_build(self):
        conda_dir = self.cwd / ".snakemake/conda"
        assert not any(conda_dir.iterdir())
