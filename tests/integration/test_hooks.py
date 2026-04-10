import pytest
from helpers import (
    ShowyourworkRepositoryActions,
    TemporaryShowyourworkRepository,
)

from showyourwork.config import edit_yaml


class TestHooks(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """
    Test the showyourwork patches to snakemake's onstart and onsuccess hooks.
    """

    hooks = ["onstart", "onsuccess"]

    def add_hook(self, name):
        with open(self.cwd / "Snakefile", "a") as f:
            print("\n", file=f)
            print(
                "\n".join(
                    [
                        f"{name}:",
                        f"    with open('{name}.txt', 'w') as f:",
                        f"        f.write('{name}\\n')",
                    ]
                ),
                file=f,
            )
            print("\n", file=f)

    def add_hooks(self):
        for hook in self.hooks:
            self.add_hook(hook)

    def customize(self):

        # Add the pipeline script
        self.add_pipeline_script()

        # Add the Snakefile rule to generate the dataset
        self.add_pipeline_rule()
        self.add_hooks()

        # Add the script to generate the figure
        self.add_figure_script(load_data=True)

        # Make the dataset a dependency of the figure
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            config["dependencies"] = {
                "src/scripts/test_figure.py": ["src/data/test_data.npz"]
            }
            config["run_cache_rules_on_ci"] = True

        # Add the figure environment to the tex file
        self.add_figure_environment()

    def check_build(self):
        for hook in self.hooks:
            hook_file = self.cwd / f"{hook}.txt"
            assert hook_file.is_file()
            assert hook_file.read_text() == f"{hook}\n"


class TestHooksError(TestHooks):
    """
    Test the showyourwork patches to snakemake's onstart and onerror hooks.
    """

    hooks = ["onstart", "onerror"]
    local_build_only = True

    def customize(self):
        super().customize()
        with open(self.cwd / "src" / "scripts" / "test_data.py", "a") as f:
            print("\nraise RuntimeError('there was an error')\n", file=f)

    def build_local(self, env=None):
        with pytest.raises(Exception, match="there was an error"):
            super().build_local(env=env)
