import pytest
from helpers import (
    ShowyourworkRepositoryActions,
    TemporaryShowyourworkRepository,
)

pytestmark = pytest.mark.remote

# Add code to the workflow to explicitly install a TeX distribution
workflow_replace = "      - name: Build the article PDF"
workflow_replace_with = """      - name: Install TinyTex for matplotlib LaTeX rendering
        id: tinytex
        shell: bash -l {0}
        run: |
          wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh" | sh
          sudo ~/bin/tlmgr install type1cm cm-super

      - name: Build the article PDF
"""


# A figure script that requires LaTeX
figure_script = r"""
# Ensure ~/bin is in the PATH so matplotlib can find latex
import os
from pathlib import Path
os.environ["PATH"] += os.pathsep + str(Path.home() / "bin")

import matplotlib.pyplot as plt
import paths
plt.rcParams.update({"text.usetex": True})
fig = plt.figure(figsize=(7, 6))
plt.plot([0, 1], [0, 1])
plt.xlabel(r'$\alpha \beta \frac{1}{2} \exp{-x^2}$')
fig.savefig(paths.figures / 'test_figure.pdf', bbox_inches='tight')
"""


class TestMatplotlibLatex(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """Test rendering TeX in matplotlib."""

    def customize(self):
        """Add a figure that requires LaTeX to build."""
        # Edit the workflow
        with open(self.cwd / ".github" / "workflows" / "build.yml", "r") as f:
            workflow = f.read().replace(
                workflow_replace, workflow_replace_with
            )
        with open(self.cwd / ".github" / "workflows" / "build.yml", "w") as f:
            print(workflow, file=f)

        # Create the figure script
        with open(self.cwd / "src" / "scripts" / "test_figure.py", "w") as f:
            print(figure_script, file=f)

        # Add the figure environment to the tex file
        self.add_figure_environment()

    def build_local(self):
        """Disable local build, since this is explicitly a CI test."""
        pass
