from helpers import (
    ShowyourworkRepositoryActions,
    TemporaryShowyourworkRepository,
)


class TestDuplicateGraphics(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """Test graphics included in multiple figure environments."""

    # No need to test this on CI
    local_build_only = True

    def customize(self):
        self.add_figure_script()
        self.add_figure_environment(label="fig:fig1")
        self.add_figure_environment(label="fig:fig2")
