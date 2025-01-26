from helpers import TemporaryShowyourworkRepository
from yaml import safe_dump, safe_load

from showyourwork.subproc import get_stdout

new_directory_script = r"""

import paths

newdir = paths.data / "newdir"

newdir.mkdir(exist_ok=True)

with open(newdir / "newfile.txt", "w") as f:
    f.write("This is a new file")
"""


class TestCleanForce(TemporaryShowyourworkRepository):
    """Test that when calling clean with the -f option a subdirectory is removed."""

    local_build_only = True

    def customize(self):
        """Create and edit all the necessary files for the workflow."""
        # Create the script
        print(f"[{self.repo}] Creating new stuff")
        with open(self.cwd / "src" / "scripts" / "new_data.py", "w") as f:
            print(new_directory_script, file=f)

        # Make the manuscript depend on this new file
        with open(self.cwd / "showyourwork.yml") as config:
            original_cfg = safe_load(config)
        with open(self.cwd / "showyourwork.yml", mode="w") as config:
            original_cfg["dependencies"] = {
                "src/tex/ms.tex": ["src/data/newdir/newfile.txt"]
            }
            safe_dump(original_cfg, config)

        # Add a Snakemake rule to run the script
        with open(self.cwd / "Snakefile") as f:
            contents = f.read()
        with open(self.cwd / "Snakefile", "w") as f:
            print(contents, file=f)
            print("\n", file=f)
            print(
                "\n".join(
                    [
                        "rule make_new_data:",
                        "    output:",
                        "        'src/data/newdir/newfile.txt'",
                        "    script:",
                        "        'src/scripts/new_data.py'",
                    ]
                ),
                file=f,
            )

    def check_build(self):
        """Check that the new directory and file were created."""
        assert (self.cwd / "src" / "data" / "newdir").exists()
        assert (self.cwd / "src" / "data" / "newdir" / "newfile.txt").exists()

        # Run the clean command with the -f option
        get_stdout("showyourwork clean -f", shell=True, cwd=self.cwd)

        # Check that the new directory was removed
        assert not (self.cwd / "src" / "data" / "newdir").exists()

class TestCleanDeep(TemporaryShowyourworkRepository):
    """Test that when calling clean with the -d option a subdirectory is removed."""

    local_build_only = True

    def check_build(self):
        # Run the clean command with the -f option
        get_stdout("showyourwork clean -d", shell=True, cwd=self.cwd)

        # Check that the new directory was removed
        assert not (self.cwd / ".snakemake").exists()
        assert not (self.cwd / ".showyourwork").exists()
