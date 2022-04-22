from helpers import TemporaryShowyourworkRepository, edit_yaml
from showyourwork import overleaf
from showyourwork.subproc import get_stdout
from tempfile import NamedTemporaryFile


class TestDefault(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo."""

    pass


class TestZenodo(TemporaryShowyourworkRepository):
    """Test a repo that downloads data from Zenodo."""

    pass


class TestOverleaf(TemporaryShowyourworkRepository):
    """Test a repo that integrates with an Overleaf project."""

    overleaf_id = "6262c032aae5421d6d945acf"

    def startup(self):
        """Wipe the Overleaf remote to start fresh."""
        overleaf.wipe_remote(self.overleaf_id)

    def customize(self):
        """Customize the build to test all Overleaf functionality."""
        # Add Overleaf settings to the config
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["overleaf"]["push"] = ["src/tex/figures"]
            config["overleaf"]["pull"] = ["src/tex/ms.tex"]

        # Mock changes on the Overleaf remote by modifying the local
        # manuscript (adding a figure environment),
        # pushing it to Overleaf, and restoring the original local version
        ms = self.cwd / "src" / "tex" / "ms.tex"
        with open(ms, "r") as f:
            ms_orig = f.read()
        with open(ms, "w") as f:
            ms_new = ms_orig.replace(
                r"\end{document}",
                "\n".join(
                    [
                        r"\begin{figure}[ht!]",
                        r"\script{random_numbers.py}",
                        r"\begin{centering}",
                        r"\includegraphics[width=\linewidth]{figures/random_numbers.pdf}",
                        r"\caption{Random numbers.}",
                        r"\label{fig:random_numbers}",
                        r"\end{centering}",
                        r"\end{figure}",
                        "",
                        r"\end{document}",
                    ]
                ),
            )
            print(ms_new, file=f)
        overleaf.push_files([ms], self.overleaf_id, path=self.cwd)
        with open(ms, "w") as f:
            print(ms_orig, file=f)

    def check_build(self):
        """Check that our figure was built correctly.

        This will only work if we successfully synced the Overleaf manuscript
        to the local repo, as that contains the figure environment defining
        the script that produces `random_numbers.pdf`.

        Here we also check that the generated figure was pushed to the
        Overleaf project successfully.
        """
        # The programmatically-generated figure
        figure = self.cwd / "src" / "tex" / "figures" / "random_numbers.pdf"

        # Check that figure is present locally
        assert (figure).exists()

        # Check that figure is present on the remote
        overleaf.pull_files(
            [figure], self.overleaf_id, path=self.cwd, error_if_missing=True
        )


class TestZenodoCache(TemporaryShowyourworkRepository):
    """Test the Zenodo caching feature."""

    zenodo_cache = True

    def customize(self):
        """Add the necessary dependencies to the config file."""
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/plot_dataset.py": "src/data/dataset.npz"
            }


class TestZenodoDirCache(TemporaryShowyourworkRepository):
    """Test the Zenodo caching feature for entire directories."""

    zenodo_cache = True

    def customize(self):
        """Add the necessary dependencies to the config file."""
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/plot_dataset.py": "src/data/dataset"
            }