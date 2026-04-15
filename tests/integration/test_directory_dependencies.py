"""Tests for directory dependencies feature.

These tests verify that dependencies can be specified as entire directories,
with all files within the directory treated as dependencies of a script.
"""

from helpers import (
    ShowyourworkRepositoryActions,
    TemporaryShowyourworkRepository,
)

from showyourwork.config import edit_yaml


class TestDirectoryDependency(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """
    Test that dependencies can be specified as directories.

    Covers three use cases:
    - directory dependency with nested subdirectories
    - multiple directory dependencies
    - mixed individual file and directory dependencies

    Also verifies that modifying a file inside a dependency directory
    triggers a rebuild of the figure.
    """

    # No need to test this on CI
    local_build_only = True

    def _create_helper_scripts(self):
        """Create helper scripts in utils/ and config/ directories."""
        utils_dir = self.cwd / "src" / "scripts" / "utils"
        utils_dir.mkdir(parents=True, exist_ok=True)

        (utils_dir / "helper1.py").write_text(
            "def helper1_function():\n    return 'Hello from helper1'\n"
        )
        (utils_dir / "helper2.py").write_text(
            "def helper2_function():\n    return 'World from helper2'\n"
        )

        nested_dir = utils_dir / "nested"
        nested_dir.mkdir(parents=True, exist_ok=True)
        (nested_dir / "helper3.py").write_text(
            "def helper3_function():\n    return 'Deep helper from helper3'\n"
        )

        config_dir = self.cwd / "src" / "scripts" / "config"
        config_dir.mkdir(parents=True, exist_ok=True)
        (config_dir / "settings.py").write_text("PLOT_SIZE = (7, 6)\n")

        (self.cwd / "src" / "scripts" / "standalone_helper.py").write_text(
            "def standalone():\n    return 'standalone'\n"
        )

    def _create_figure_script(self, seed=0):
        """Create the main figure script that imports from helpers."""
        with open(self.cwd / "src" / "scripts" / "test_figure.py", "w") as f:
            f.write(
                "import sys\n"
                "from pathlib import Path\n"
                "sys.path.insert(0, str(Path.cwd() / 'src' / 'scripts'))\n"
                "from utils.helper1 import helper1_function\n"
                "from utils.helper2 import helper2_function\n"
                "from utils.nested.helper3 import helper3_function\n"
                "from config.settings import PLOT_SIZE\n"
                "from standalone_helper import standalone\n"
                "import matplotlib.pyplot as plt\n"
                "import numpy as np\n"
                "import paths\n"
                f"np.random.seed({seed})\n"
                "msg1 = helper1_function()\n"
                "msg2 = helper2_function()\n"
                "msg3 = helper3_function()\n"
                "msg4 = standalone()\n"
                "fig = plt.figure(figsize=PLOT_SIZE)\n"
                "plt.plot(np.random.randn(10))\n"
                "plt.title(f'{msg1} {msg2} {msg3} {msg4}')\n"
                "fig.savefig(paths.figures / 'test_figure.pdf', "
                "bbox_inches='tight', dpi=300)\n"
            )

    def _configure_dependencies(self):
        """Set up mixed file and directory dependencies in showyourwork.yml."""
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/test_figure.py": [
                    "src/scripts/utils/",  # Directory dependency (with nested subdir)
                    "src/scripts/config/",  # Second directory dependency
                    "src/scripts/standalone_helper.py",  # Individual file dependency
                ],
            }

    def customize(self):
        """Create figure script, helper scripts, and dependencies config."""
        self._create_helper_scripts()
        self._create_figure_script(seed=0)
        self.add_figure_environment()
        self._configure_dependencies()

    def check_build(self):
        """Verify build and that modifying a dependency triggers a rebuild."""
        figure_path = self.cwd / "src" / "tex" / "figures" / "test_figure.pdf"
        assert (
            figure_path.exists()
        ), "Figure was not generated with directory dependencies"

        # Record the modification time of the generated figure
        mtime_before = figure_path.stat().st_mtime

        # Modify a helper inside a dependency directory
        (self.cwd / "src" / "scripts" / "utils" / "helper1.py").write_text(
            "def helper1_function():\n    return 'Modified helper1'\n"
        )

        # Commit and rebuild
        self.git_commit()
        self.build_local()

        # The figure should have been regenerated
        mtime_after = figure_path.stat().st_mtime
        assert (
            mtime_after > mtime_before
        ), "Figure was not regenerated after modifying a dependency inside a directory"
