from helpers import TemporaryShowyourworkRepository
from helpers.temp_repo import SANDBOX


class TestSpaceInPath(TemporaryShowyourworkRepository):
    """Test for issues when local path has a space."""

    local_build_only = True
    root_path = SANDBOX / "space in path"
