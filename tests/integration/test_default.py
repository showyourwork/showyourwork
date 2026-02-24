from helpers import TemporaryShowyourworkRepository


class TestDefault(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo."""

    pass


class TestDefaultDry(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo with a dry-run"""

    dry_run = True
