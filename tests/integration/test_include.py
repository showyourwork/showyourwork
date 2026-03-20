import shutil

from helpers import ShowyourworkRepositoryActions, TemporaryShowyourworkRepository

from showyourwork.config import edit_yaml


class TestInclude(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """Test a workflow with rules split accross files"""

    def customize(self):
        """Add all files for the build"""
        # Add the pipeline python script
        self.add_pipeline_script()

        # Create a snakefile
        self.add_pipeline_rule()

        # Move the snakefile to a rule file
        snakefile = self.cwd / "Snakefile"
        rulefile = self.cwd / "src/rules/test_data.smk"
        rulefile.parent.mkdir()
        shutil.move(snakefile, rulefile)

        # Update the script path in the rule file
        rule_txt = rulefile.read_text()
        rule_txt = rule_txt.replace(
            "'src/scripts/test_data.py'",
            "'../../src/scripts/test_data.py'",
        )
        rulefile.write_text(rule_txt)

        # Snakefile now only includes the rule
        with open(snakefile, "w") as f:
            print("include: 'src/rules/test_data.smk'", file=f)

        # Add the script to generate the figure
        self.add_figure_script(load_data=True)

        # Make the dataset a dependency of the figure
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/test_figure.py": "src/data/test_data.npz"
            }
            config["run_cache_rules_on_ci"] = True

        # Add the figure environment to the tex file
        self.add_figure_environment()
