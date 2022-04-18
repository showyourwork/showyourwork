from helpers import TemporaryShowyourworkRepository
from pathlib import Path

'''
class TestDefault(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo."""

    pass


class TestDefaultWithZenodoCaching(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo with Zenodo caching."""

    def create_local(self):
        super().create_local(zenodo_cache=True)
'''


class TestZenodo(TemporaryShowyourworkRepository):
    """Test a repo that downloads data from Zenodo."""

    pass