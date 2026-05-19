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
            micromamba_path = os.getenv("MAMBA_EXE", "micromamba")
            conda_prefix = os.getenv("CONDA_PREFIX")

            if conda_prefix:
                install_target = f'-p "{conda_prefix}"'
            else:
                # Fall back to micromamba's active/default target if no prefix is set.
                install_target = ""

            get_stdout(
                f"{micromamba_path} install -y {install_target} tectonic=0.14.1",
                cwd=self.cwd,
                shell=True,
            )
            get_stdout(
                (
                    f"{micromamba_path} run {install_target} "
                    "conda pypi install numpy matplotlib"
                ),
                cwd=self.cwd,
                shell=True,
            )

    def build_local(self):
        super().build_local(args=["--no-conda"])

    def check_build(self):
        conda_dir = self.cwd / ".snakemake/conda"
        assert not any(conda_dir.iterdir())
