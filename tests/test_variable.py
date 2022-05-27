from helpers import TemporaryShowyourworkRepository
from showyourwork.config import edit_yaml


# A script that computes the age of the universe
variable_script = r"""
import paths
import numpy as np

# Compute the age of the universe
np.random.seed(42)
age = np.random.normal(14.0, 1.0)

# Write it to disk
with open(paths.output / "age_of_universe.txt", "w") as f:
    f.write(f"{age:.3f}")
"""


class TestLatexVariable(TemporaryShowyourworkRepository):
    """Test a workflow with dynamic quantities imported into the tex file."""

    # No need to test this on CI
    local_build_only = True

    def customize(self):
        """Create and edit all the necessary files for the workflow."""
        # Create the script
        with open(
            self.cwd / "src" / "scripts" / "age_of_universe.py", "w"
        ) as f:
            print(variable_script, file=f)

        # Import the variable into the tex file
        ms = self.cwd / "src" / "tex" / "ms.tex"
        with open(ms, "r") as f:
            ms_orig = f.read()
        with open(ms, "w") as f:
            ms_new = ms_orig.replace(
                r"\end{document}",
                r"Based on a detailed analysis of Planck observations of the cosmic "
                r"microwave background, we have determined the age of the universe "
                r"to be \variable{output/age_of_universe.txt} Gyr."
                "\n"
                r"\end{document}",
            )
            print(ms_new, file=f)

        # Add a Snakemake rule to run the script
        with open(self.cwd / "Snakefile", "r") as f:
            contents = f.read()
        with open(self.cwd / "Snakefile", "w") as f:
            print(contents, file=f)
            print("\n", file=f)
            print(
                "\n".join(
                    [
                        "rule age_of_universe:",
                        "    output:",
                        "        'src/tex/output/age_of_universe.txt'",
                        "    script:",
                        "        'src/scripts/age_of_universe.py'",
                    ]
                ),
                file=f,
            )