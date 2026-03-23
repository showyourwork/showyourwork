from helpers import TemporaryShowyourworkRepository


class TestDefault(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo."""

    pass


class TestDefaultDry(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo with a dry-run"""

    def build_local(self):
        super().build_local(args=["--dry-run"])
