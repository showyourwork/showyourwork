"""Tests about handling figures and graphics."""

import shutil

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


class TestStaticFigures(TemporaryShowyourworkRepository, ShowyourworkRepositoryActions):
    """Test aaddition of static figures."""

    # No need to test this on CI
    local_build_only = True

    def customize(self):
        shutil.copy(
            self.cwd / "src/tex/orcid-ID.png", self.cwd / "src/static/orcid-ID.png"
        )

        (self.cwd / "src/static/a/b").mkdir(exist_ok=True, parents=True)

        shutil.copy(
            self.cwd / "src/tex/orcid-ID.png", self.cwd / "src/static/a/orcid-ID.png"
        )

        shutil.copy(
            self.cwd / "src/tex/orcid-ID.png",
            self.cwd / "src/static/a/b/orcid-ID.png",
        )

        self.add_figure_environment(
            figure_path="orcid-ID.png",
            label="fig:orcid-id",
            add_script=False,
        )
        self.add_figure_environment(
            figure_path="a/orcid-ID.png",
            label="fig:orcid-id-2",
            add_script=False,
        )
        self.add_figure_environment(
            figure_path="a/b/orcid-ID.png",
            label="fig:orcid-id-3",
            add_script=False,
        )

    def check_build(self):
        """Check that static images have been copied to figures respecting directory structure."""
        assert (self.cwd / "src/tex/figures/orcid-ID.png").exists()
        assert (self.cwd / "src/tex/figures/a/orcid-ID.png").exists()
        assert (self.cwd / "src/tex/figures/a/b/orcid-ID.png").exists()
