from helpers import TemporaryShowyourworkRepository
from showyourwork import overleaf, exceptions
from showyourwork.config import edit_yaml
from showyourwork.logging import get_logger
from showyourwork.subproc import get_stdout
from tempfile import NamedTemporaryFile
import time


class TestDefault(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo."""

    pass


class TestZenodo(TemporaryShowyourworkRepository):
    """Test a repo that downloads data from Zenodo."""

    def customize(self):
        """Add the necessary dependencies to the config file."""
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/TOI640.py": [
                    "src/data/TOI640/planet.json",
                    "src/data/TOI640/S06.txt",
                    "src/data/TOI640/S06.json",
                    "src/data/TOI640/S07.txt",
                    "src/data/TOI640/S07.json",
                ]
            }
            config["datasets"] = {
                "10.5281/zenodo.6515864": {
                    "destination": "src/data/TOI640",
                    "contents": {
                        "README.md": None,
                        "TOI640b.json": "src/data/TOI640/planet.json",
                        "images.tar.gz": {
                            "README.md": None,
                            "S06": {"image.json": "src/data/TOI640/S06.json"},
                            "S07": {"image.json": "src/data/TOI640/S07.json"},
                        },
                        "lightcurves.zip": {
                            "lightcurves": {
                                "README.md": None,
                                "S06": {"lc.txt": "src/data/TOI640/S06.txt"},
                                "S07": {"lc.txt": "src/data/TOI640/S07.txt"},
                            }
                        },
                    },
                }
            }


class TestOverleaf(TemporaryShowyourworkRepository):
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
        for n in range(self.auth_retries):
            try:
                overleaf.push_files([ms], self.overleaf_id, path=self.cwd)
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
        `random_numbers.pdf`.

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
        figure = self.cwd / "src" / "tex" / "figures" / "random_numbers.pdf"
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


class TestSandboxCache(TemporaryShowyourworkRepository):
    """Test the Zenodo caching feature."""

    # Enable Sandbox caching
    sandbox_cache = True

    def customize(self):
        """Add the necessary dependencies to the config file."""
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/plot_dataset.py": "src/data/dataset.npz"
            }


class TestSandboxDirCache(TemporaryShowyourworkRepository):
    """Test the Zenodo caching feature for entire directories."""

    # Enable Sandbox caching
    sandbox_cache = True

    def customize(self):
        """Add the necessary dependencies to the config file."""
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/plot_dataset.py": "src/data/dataset"
            }


class TestFallbackRules(TemporaryShowyourworkRepository):
    """Test the fallback rules that inform the user about missing article dependencies."""

    # No need to test this on CI
    local_build_only = True

    def add_figure_environment(self):
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

    def add_figure_script(self):
        with open(
            self.cwd / "src" / "scripts" / "random_numbers.py", "w"
        ) as f:
            print(
                "\n".join(
                    [
                        "import matplotlib.pyplot as plt",
                        "import numpy as np",
                        "import paths",
                        "random_numbers = np.random.randn(100, 10)",
                        "fig = plt.figure(figsize=(7, 6))",
                        "plt.plot(random_numbers)",
                        "fig.savefig(paths.figures / 'random_numbers.pdf', bbox_inches='tight', dpi=300)",
                    ]
                ),
                file=f,
            )

    def build_local(self):
        # Build a blank article
        super().build_local()

        # Add a figure environment to the tex file w/ a script command,
        # but don't actually create the script. The build should fail
        # with an informative message.
        self.add_figure_environment()
        try:
            super().build_local()
        except Exception as e:
            if (
                not "No rule to generate src/tex/figures/random_numbers.pdf"
                in str(e)
            ):
                raise Exception("Incorrect exception message.")
        else:
            raise Exception(
                "Expected failure, but article build step succeeded."
            )

        # Now add the script. The build should succeed
        self.add_figure_script()
        super().build_local()