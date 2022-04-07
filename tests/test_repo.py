from temp_repo import TemporaryShowyourworkRepository


class TestDefault(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo."""

    pass


class TestDefaultWithZenodoCaching(TemporaryShowyourworkRepository):
    """Test setting up and building the default repo with Zenodo caching."""

    def create_local(self):
        super().create_local(zenodo_cache=True)