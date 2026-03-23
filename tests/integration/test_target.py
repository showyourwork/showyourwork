from helpers import ShowyourworkRepositoryActions, TemporaryShowyourworkRepository


class TestTarget(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """The setting up and building the repo with target rule"""

    def customize(self):
        self.add_pipeline_script()

        self.add_pipeline_rule()

    def build_local(self):
        super().build_local(args=["generate_data"])
