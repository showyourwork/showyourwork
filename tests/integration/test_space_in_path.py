import shutil

from helpers import ShowyourworkRepositoryActions, TemporaryShowyourworkRepository
from helpers.temp_repo import SANDBOX


class TestSpaceInPath(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """Test for issues when local path has a space."""

    local_build_only = True
    root_path = SANDBOX / "space in path"

    def customize(self):
        # Static figures are copied with a shell command containing the path
        shutil.copy(
            self.cwd / "src/tex/orcid-ID.png", self.cwd / "src/static/orcid-ID.png"
        )
        self.add_figure_environment(
            figure_path="orcid-ID.png",
            label="fig:orcid-id",
            add_script=False,
        )

    def check_build(self):
        assert (self.cwd / "src/tex/figures/orcid-ID.png").exists()
