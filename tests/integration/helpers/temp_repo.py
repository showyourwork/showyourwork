import asyncio
import inspect
import logging
import os
import re
import shutil
from pathlib import Path

import pytest
import yaml

import showyourwork
from showyourwork import gitapi
from showyourwork.config import render_config
from showyourwork.git import get_repo_sha
from showyourwork.logging import get_logger
from showyourwork.subproc import get_stdout
from showyourwork.zenodo import Zenodo

try:
    from yaml import CDumper as Dumper
    from yaml import CLoader as Loader
except ImportError:
    # If LibYAML not installed
    from yaml import Dumper, Loader


# Folder where we'll create our temporary repos
SANDBOX = Path(__file__).absolute().parents[1] / "sandbox"
RESOURCES = Path(__file__).absolute().parents[1] / "resources"


class TemporaryShowyourworkRepository:
    """
    Base class for showyourwork tests.

    """

    # Local debug mode
    debug = os.getenv("DEBUG", "False").lower() in ("true", "1", "t")

    # Editable class settings
    cache = False
    overleaf_id = None
    action_wait = 240
    action_max_tries = 10
    action_interval = 60
    use_local_showyourwork = debug
    showyourwork_version = None
    local_build_only = debug
    require_local_build = False
    delete_remote_on_success = False
    clear_actions_cache_on_start = True

    # Internal
    _sandbox_concept_doi = None

    @property
    def repo(self):
        # Name the repo after the subclass (in hyphenated-snake-case)
        return (
            re.sub(r"(?<!^)(?=[A-Z])", "_", self.__class__.__name__)
            .lower()
            .replace("_", "-")
        )

    @property
    def cwd(self):
        return SANDBOX / self.repo

    def startup(self):
        """Subclass me to run things at startup."""
        pass

    def teardown(self):
        """Subclass me to run things at teardown."""
        pass

    def customize(self):
        """Subclass me to customize the repo."""
        pass

    def check_build(self):
        """Subclass me to run post-build local checks."""
        pass

    def create_local(self, use_local=False):
        """Create the repo locally."""
        # Delete any local repos
        self.delete_local()

        # Parse options
        command = "showyourwork setup"
        if use_local or self.showyourwork_version is None:
            if use_local or self.use_local_showyourwork:
                version = str(Path(showyourwork.__file__).parents[1])
            else:
                version = get_repo_sha()
        else:
            version = self.showyourwork_version
        options = f"--quiet --version={version} "
        if self.cache:
            # Enable zenodo caching
            options += "--cache"
        if self.overleaf_id:
            # Enable overleaf syncing
            options += f"--overleaf={self.overleaf_id}"

        # Provide git name & email
        if os.getenv("CI", "false") == "true":
            get_stdout(
                "git config --global user.name 'gh-actions'",
                cwd=SANDBOX,
                shell=True,
            )
            get_stdout(
                "git config --global user.email 'gh-actions'",
                cwd=SANDBOX,
                shell=True,
            )

        # Create a new one
        print(
            f"[{self.repo}] Creating local repo `tests/sandbox/{self.repo}`..."
        )
        get_stdout(
            f"{command} {options} showyourwork/{self.repo}",
            cwd=SANDBOX,
            shell=True,
        )

        # Get the Zenodo Sandbox cache concept doi for the main branch (if any)
        self._sandbox_concept_doi = (
            render_config(cwd=self.cwd)
            .get("cache", {})
            .get("main", {})
            .get("sandbox")
        )

        # Tweak the README
        message = "\n".join(
            [line.strip() for line in self.__doc__.split("\n")]
        )
        file = Path(inspect.getfile(self.__class__)).name
        readme = (
            "*This is an automatically generated test for "
            "[showyourwork](https://github.com/showyourwork/showyourwork) "
            "generated from the file "
            f"[{file}](https://github.com/showyourwork/showyourwork/blob/main/tests/integration/{file}).*"
            "\n\n" + message
        )
        with open(self.cwd / "README.md", "r") as f:
            contents = f.read()
        contents = contents.replace(
            "An open source scientific article created using the "
            "[showyourwork](https://github.com/showyourwork/showyourwork) "
            "workflow.",
            readme,
        )
        with open(self.cwd / "README.md", "w") as f:
            f.write(contents)

        # Update the version in the config file to point to the correct fork if required
        fork = os.environ.get("SHOWYOURWORK_FORK", "").strip()
        if fork:
            ref = (
                os.environ.get("SHOWYOURWORK_REF", "").strip()
                or get_repo_sha()
            )
            with open(self.cwd / "showyourwork.yml", "r") as f:
                contents = yaml.load(f, Loader=Loader)
            contents["version"] = {"fork": fork, "ref": ref}
            with open(self.cwd / "showyourwork.yml", "w") as f:
                yaml.dump(contents, f, Dumper=Dumper)

    def create_remote(self):
        """Create the repo on GitHub if needed."""
        print(
            f"[{self.repo}] Setting up remote GitHub repo "
            f"`showyourwork/{self.repo}`..."
        )
        gitapi.create_repo(
            self.repo,
            org="showyourwork",
            description="Temporary test repository for showyourwork",
            private=False,
        )

    def clear_remote_cache(self):
        """Clear the remote GitHub Actions cache."""
        print(f"[{self.repo}] Clearing the Actions cache...")
        gitapi.clear_cache(self.repo, org="showyourwork")

    def git_commit(self):
        """Add and commit all files in the local repo."""
        get_stdout("git add .", shell=True, cwd=self.cwd)
        get_stdout(
            "git diff-index --quiet HEAD || "
            "git -c user.name='gh-actions' -c user.email='gh-actions' "
            "commit -q -m 'auto commit from showyourwork tests'",
            shell=True,
            cwd=self.cwd,
        )

    def build_local(self, pre=""):
        """Run showyourwork locally to build the article."""
        print(f"[{self.repo}] Building the article locally...")

        def callback(code, stdout, stderr):
            if code != 0:
                raise Exception(stdout + "\n" + stderr)

        get_stdout(
            f"{pre} CI=false showyourwork build",
            shell=True,
            cwd=self.cwd,
            callback=callback,
        )

    @pytest.mark.asyncio_cooperative
    async def run_github_action(self):
        """
        Push to the remote and asynchronously wait for the workflow on
        GitHub Actions to finish.
        On success, returns, otherwise raises an Exception.

        """
        print(f"[{self.repo}] Pushing to `showyourwork/{self.repo}`...")
        get_stdout(
            "git push --force https://x-access-token:"
            f"{gitapi.get_access_token()}"
            f"@github.com/showyourwork/{self.repo} main",
            shell=True,
            cwd=self.cwd,
            secrets=[gitapi.get_access_token()],
        )
        head_sha = get_stdout(
            "git rev-parse HEAD", shell=True, cwd=self.cwd
        ).replace("\n", "")
        print(
            f"[{self.repo}] Waiting {self.action_wait} seconds for workflow "
            f"to finish (1/{self.action_max_tries})..."
        )
        await asyncio.sleep(self.action_wait)
        status = "unknown"
        for n in range(self.action_max_tries):
            (status, conclusion, url,) = gitapi.get_workflow_run_status(
                self.repo,
                org="showyourwork",
                q={"event": "push", "head_sha": head_sha},
            )
            if status == "completed":
                if conclusion == "success":
                    print(f"[{self.repo}] Workflow completed successfully.")
                    return
                else:
                    raise Exception(
                        "[{self.repo}] GitHub Actions workflow terminated "
                        f"with status {conclusion}.\n"
                        f"For details, see {url}"
                    )
            elif n < self.action_max_tries - 1:
                print(
                    f"[{self.repo}] Waiting {self.action_interval} seconds for "
                    f"workflow to finish ({n+2}/{self.action_max_tries})..."
                )
                await asyncio.sleep(self.action_interval)
        else:
            raise Exception(
                "[{self.repo}] GitHub Actions workflow timed out.\n"
                f"For details, see {url}"
            )

    def delete_zenodo(self):
        """Delete the Zenodo deposit associated with the temp repo."""
        if self._sandbox_concept_doi:
            print(
                f"[{self.repo}] Deleting Zenodo Sandbox deposit "
                f"with concept DOI {self._sandbox_concept_doi}..."
            )
            Zenodo(self._sandbox_concept_doi).delete()

    def delete_remote(self):
        """Delete the remote repo."""
        print(
            f"[{self.repo}] Deleting GitHub repo "
            f"`showyourwork/{self.repo}`..."
        )
        gitapi.delete_repo(self.repo, org="showyourwork")

    def delete_local(self):
        """Delete the local repo."""
        if (self.cwd).exists():
            print(
                f"[{self.repo}] Deleting local repo `tests/sandbox/{self.repo}`..."
            )
            shutil.rmtree(self.cwd)

    def disable_logging(self):
        """Disable showyourwork screen output."""
        for handler in get_logger().handlers:
            if isinstance(handler, logging.StreamHandler):
                handler.setLevel(logging.ERROR)

    def test_local(self):
        """
        Test functionality by creating a new repo, customizing it,
        pushing it to GitHub, and awaiting the article build action to
        complete.

        """
        try:

            # Disable screen logging info from showyourwork
            self.disable_logging()

            # Always run this first
            self.startup()

            # Set up the repo
            self.create_local(use_local=True)

            # Customize the repo
            self.customize()

            # Commit changes
            self.git_commit()

            # Build the article locally
            self.build_local()

            # Run local checks
            self.check_build()

            # Delete local repo (only on success)
            self.delete_local()

        finally:

            # Always delete the Zenodo deposit (if created)
            self.delete_zenodo()

            # Always run this last
            self.teardown()

    @pytest.mark.remote
    @pytest.mark.asyncio_cooperative
    async def test_remote(self):
        """
        Test functionality by creating a new repo, customizing it,
        pushing it to GitHub, and awaiting the article build action to
        complete.

        """
        if self.local_build_only:
            return

        try:

            # Disable screen logging info from showyourwork
            self.disable_logging()

            # Always run this first
            self.startup()

            # Create the repo on GitHub
            self.create_remote()
            if self.clear_actions_cache_on_start:
                self.clear_remote_cache()

            # Set up the repo
            self.create_local()

            # Customize the repo
            self.customize()

            # Commit changes
            self.git_commit()

            # Some tests require a local build first
            if self.require_local_build:
                self.build_local()
                self.check_build()

            # Push to GitHub to trigger the Actions workflow
            # and wait for the result
            await self.run_github_action()

            # Delete remote repo (only on success)
            if self.delete_remote_on_success:
                self.delete_remote()

            # Delete local repo (only on success)
            self.delete_local()

        finally:

            # Always delete the Zenodo deposit (if created)
            self.delete_zenodo()

            # Always run this last
            self.teardown()


class ShowyourworkRepositoryActions:
    def add_figure_environment(self, add_script=True, label="fig:test_figure"):
        """Adds a figure environment to the TeX file that includes `test_figure.pdf`."""
        ms = self.cwd / "src" / "tex" / "ms.tex"
        with open(ms, "r") as f:
            ms_orig = f.read()
        with open(ms, "w") as f:
            ms_new = ms_orig.replace(
                r"\end{document}",
                "\n".join(
                    [
                        r"\begin{figure}[ht!]",
                        r"\script{test_figure.py}" if add_script else "",
                        r"\begin{centering}",
                        r"\includegraphics[width=\linewidth]{figures/test_figure.pdf}",
                        r"\caption{A test figure.}",
                        r"\label{" + label + "}",
                        r"\end{centering}",
                        r"\end{figure}",
                        "",
                        r"\end{document}",
                    ]
                ),
            )
            print(ms_new, file=f)

    def add_figure_script(self, load_data=False, batch=False):
        """Creates a figure script `test_figure.py` that generates `test_figure.pdf`."""
        if load_data:
            if batch:
                get_data = "data = np.array([np.load(paths.data / 'test_data' / f'test_data{n:02d}.npz')['data'] for n in range(50)])"
            else:
                get_data = "data = np.load(paths.data / 'test_data.npz')"
        else:
            get_data = "data = np.random.randn(100, 10)"
        with open(self.cwd / "src" / "scripts" / "test_figure.py", "w") as f:
            print(
                "\n".join(
                    [
                        "import matplotlib.pyplot as plt",
                        "import numpy as np",
                        "import paths",
                        "np.random.seed(0)",
                        get_data,
                        "fig = plt.figure(figsize=(7, 6))",
                        "plt.plot(data)",
                        "fig.savefig(paths.figures / 'test_figure.pdf', bbox_inches='tight', dpi=300)",
                    ]
                ),
                file=f,
            )

    def add_pipeline_script(self, batch=False, seed=0):
        """Creates a pipeline script `test_data.py` that generates test data."""
        if batch:
            gen_data = "\n".join(
                [
                    "(paths.data / 'test_data').mkdir(exist_ok=True)",
                    "for n in range(50):",
                    "    np.random.seed(n)",
                    "    data = np.random.randn(100)",
                    "    np.savez(paths.data / 'test_data' / f'test_data{n:02d}.npz', data=data)",
                ]
            )
        else:
            gen_data = "\n".join(
                [
                    "data = np.random.randn(100, 10)",
                    "np.savez(paths.data / 'test_data.npz', data=data)",
                ]
            )
        with open(self.cwd / "src" / "scripts" / "test_data.py", "w") as f:
            print(
                "\n".join(
                    [
                        "import numpy as np",
                        "import paths",
                        "import os",
                        "if os.getenv('CI', 'false') == 'true' or os.getenv('SYW_NO_RUN', 'false') == 'true':",
                        "    raise Exception('Output should have been downloaded from Zenodo.')",
                        f"np.random.seed({seed})",
                        gen_data,
                    ]
                ),
                file=f,
            )

    def add_pipeline_rule(self, batch=False):
        """Adds a Snakemake rule to generate test data from `test_data.py`."""
        with open(self.cwd / "Snakefile", "r") as f:
            contents = f.read()
        with open(self.cwd / "Snakefile", "w") as f:
            print(contents, file=f)
            print("\n", file=f)
            print(
                "\n".join(
                    [
                        "rule generate_data:",
                        "    output:",
                        "        directory('src/data/test_data')"
                        if batch
                        else "        'src/data/test_data.npz'",
                        "    cache:",
                        "        True",
                        "    script:",
                        "        'src/scripts/test_data.py'",
                    ]
                ),
                file=f,
            )
