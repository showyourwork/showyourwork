"""Tests for directory dependencies feature.

These tests verify that dependencies can be specified as entire directories,
with all files within the directory treated as dependencies of a script.

The tests are designed to be cross-platform compatible (Windows, macOS, Linux)
by using pathlib for filesystem operations and forward slashes (/) in YAML
configuration paths.
"""

import yaml
from helpers import (
    ShowyourworkRepositoryActions,
    TemporaryShowyourworkRepository,
)


class TestDirectoryDependency(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """
    Test that a dependency can be specified as a directory.

    When a directory is specified as a dependency in the showyourwork.yml,
    all files within that directory should be treated as dependencies of the
    dict key. This allows users to specify the entire
    utils directory rather than listing each helper script individually.
    """

    # No need to test this on CI
    local_build_only = True

    def customize(self):
        """Create figure script, helper scripts, and dependencies config."""

        # Create helper scripts in a utils subdirectory
        utils_dir = self.cwd / "src" / "scripts" / "utils"
        utils_dir.mkdir(parents=True, exist_ok=True)

        # Create first helper script
        (utils_dir / "helper1.py").write_text(
            "def helper1_function():\n    return 'Hello from helper1'\n"
        )

        # Create second helper script
        (utils_dir / "helper2.py").write_text(
            "def helper2_function():\n    return 'World from helper2'\n"
        )

        # Create a nested helper script
        nested_dir = utils_dir / "nested"
        nested_dir.mkdir(parents=True, exist_ok=True)
        (nested_dir / "helper3.py").write_text(
            "def helper3_function():\n    return 'Deep helper from helper3'\n"
        )

        # Create the main figure script that imports from utils
        with open(self.cwd / "src" / "scripts" / "test_figure.py", "w") as f:
            f.write(
                "import sys\n"
                "from pathlib import Path\n"
                "sys.path.insert(0, str(Path.cwd() / 'src' / 'scripts'))\n"
                "from utils.helper1 import helper1_function\n"
                "from utils.helper2 import helper2_function\n"
                "from utils.nested.helper3 import helper3_function\n"
                "import matplotlib.pyplot as plt\n"
                "import numpy as np\n"
                "import paths\n"
                "np.random.seed(0)\n"
                "# Use the helper functions\n"
                "msg1 = helper1_function()\n"
                "msg2 = helper2_function()\n"
                "msg3 = helper3_function()\n"
                "# Create a simple figure\n"
                "fig = plt.figure(figsize=(7, 6))\n"
                "plt.plot([1, 2, 3, 4], [1, 4, 2, 3])\n"
                "plt.title(f'{msg1} {msg2} {msg3}')\n"
                "fig.savefig(paths.figures / 'test_figure.pdf', "
                "bbox_inches='tight', dpi=300)\n"
            )

        # Add figure environment to the TeX file
        self.add_figure_environment()

        # Add dependency on the entire utils directory
        # This is the key feature being tested - specifying a directory
        # instead of individual files
        config_path = self.cwd / "showyourwork.yml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        config["dependencies"] = {
            "src/scripts/test_figure.py": [
                "src/scripts/utils/",  # Specify the entire utils directory
            ],
        }

        with open(config_path, "w") as f:
            yaml.dump(config, f, default_flow_style=False)

    def check_build(self):
        """Verify that the figure was generated successfully."""
        # The figure should have been generated
        assert (
            self.cwd / "src" / "tex" / "figures" / "test_figure.pdf"
        ).exists(), "Figure was not generated despite having directory dependency"


class TestDirectoryDependencyMultiple(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """
    Test that multiple directory dependencies work correctly.

    This verifies that a script can depend on multiple directories,
    and all files within each directory are properly tracked.
    """

    # No need to test this on CI
    local_build_only = True

    def customize(self):
        """Create multiple directories with helper scripts."""

        # Create first utils directory
        utils_dir = self.cwd / "src" / "scripts" / "utils"
        utils_dir.mkdir(parents=True, exist_ok=True)
        (utils_dir / "helper.py").write_text("def helper_func():\n    return 'data'\n")

        # Create a config directory
        config_dir = self.cwd / "src" / "scripts" / "config"
        config_dir.mkdir(parents=True, exist_ok=True)
        (config_dir / "settings.py").write_text("PLOT_SIZE = (7, 6)\n")

        # Create the main figure script
        with open(self.cwd / "src" / "scripts" / "test_figure.py", "w") as f:
            f.write(
                "import sys\n"
                "from pathlib import Path\n"
                "sys.path.insert(0, str(Path.cwd() / 'src' / 'scripts'))\n"
                "from utils.helper import helper_func\n"
                "from config.settings import PLOT_SIZE\n"
                "import matplotlib.pyplot as plt\n"
                "import numpy as np\n"
                "import paths\n"
                "np.random.seed(0)\n"
                "data = helper_func()\n"
                "fig = plt.figure(figsize=PLOT_SIZE)\n"
                "plt.plot([1, 2, 3])\n"
                "fig.savefig(paths.figures / 'test_figure.pdf', "
                "bbox_inches='tight', dpi=300)\n"
            )

        # Add figure environment
        self.add_figure_environment()

        # Add dependencies on multiple directories
        config_path = self.cwd / "showyourwork.yml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        config["dependencies"] = {
            "src/scripts/test_figure.py": [
                "src/scripts/utils/",
                "src/scripts/config/",  # Multiple directories
            ],
        }

        with open(config_path, "w") as f:
            yaml.dump(config, f, default_flow_style=False)

    def check_build(self):
        """Verify that the figure was generated successfully."""
        assert (
            self.cwd / "src" / "tex" / "figures" / "test_figure.pdf"
        ).exists(), "Figure was not generated with multiple directory dependencies"


class TestMixedFilesAndDirectories(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """
    Test that file and directory dependencies can be mixed.

    A script can depend on both individual files and entire directories
    in the same dependency list.
    """

    # No need to test this on CI
    local_build_only = True

    def customize(self):
        """Create both individual files and directories."""

        # Create a standalone helper file
        (self.cwd / "src" / "scripts").mkdir(parents=True, exist_ok=True)
        (self.cwd / "src" / "scripts" / "standalone_helper.py").write_text(
            "def standalone():\n    return 'standalone'\n"
        )

        # Create a utils directory with multiple helpers
        utils_dir = self.cwd / "src" / "scripts" / "utils"
        utils_dir.mkdir(parents=True, exist_ok=True)
        (utils_dir / "helper.py").write_text("def from_utils():\n    return 'utils'\n")

        # Create the main figure script
        with open(self.cwd / "src" / "scripts" / "test_figure.py", "w") as f:
            f.write(
                "import sys\n"
                "from pathlib import Path\n"
                "sys.path.insert(0, str(Path.cwd() / 'src' / 'scripts'))\n"
                "from standalone_helper import standalone\n"
                "from utils.helper import from_utils\n"
                "import matplotlib.pyplot as plt\n"
                "import numpy as np\n"
                "import paths\n"
                "np.random.seed(0)\n"
                "msg1 = standalone()\n"
                "msg2 = from_utils()\n"
                "fig = plt.figure(figsize=(7, 6))\n"
                "plt.plot([1, 2, 3])\n"
                "fig.savefig(paths.figures / 'test_figure.pdf', "
                "bbox_inches='tight', dpi=300)\n"
            )

        # Add figure environment
        self.add_figure_environment()

        # Add dependencies mixing files and directories
        config_path = self.cwd / "showyourwork.yml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        config["dependencies"] = {
            "src/scripts/test_figure.py": [
                "src/scripts/standalone_helper.py",  # Individual file
                "src/scripts/utils/",  # Directory
            ],
        }

        with open(config_path, "w") as f:
            yaml.dump(config, f, default_flow_style=False)

    def check_build(self):
        """Verify that the figure was generated successfully."""
        assert (
            self.cwd / "src" / "tex" / "figures" / "test_figure.pdf"
        ).exists(), (
            "Figure was not generated with mixed file and directory dependencies"
        )
