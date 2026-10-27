from helpers import TemporaryShowyourworkRepository


class TestDefault(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo."""

    def build_local(self, env=None, args=None):
        stdout = super().build_local(env=env, args=args)
        assert "Zenodo cache not yet published for this repository" not in stdout


class TestDefaultDry(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo with a dry-run"""

    def build_local(self):
        super().build_local(args=["--dry-run"])
