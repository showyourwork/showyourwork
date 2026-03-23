from helpers import ShowyourworkRepositoryActions, TemporaryShowyourworkRepository


# TODO: Add assertions that the output files exist
class TestTargetRule(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """The setting up and building the repo with target rule"""

    def customize(self):
        self.add_pipeline_script()

        self.add_pipeline_rule()

    def build_local(self):
        super().build_local(args=["generate_data"])


class TestTargetFile(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """The setting up and building the repo with target rule"""

    def customize(self):
        self.add_pipeline_script()

        self.add_pipeline_rule()

    def build_local(self):
        super().build_local(args=["src/data/test_data.npz"])


class TestTargetMulti(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """The setting up and building the repo with target rule"""

    def customize(self):
        self.add_pipeline_script()

        self.add_pipeline_rule()

        with open(self.cwd / "Snakefile", "a") as f:
            print("\n", file=f)
            print(
                "\n".join(
                    [
                        "rule touchfile:",
                        "    output:",
                        "        'src/data/afile.txt'",
                        "    shell:",
                        "        '''",
                        "        touch {output}",
                        "        '''",
                    ]
                ),
                file=f,
            )

    def build_local(self):
        super().build_local(args=["generate_data", "touchfile"])
