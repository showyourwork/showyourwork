import gzip
from pathlib import Path

from helpers import TemporaryShowyourworkRepository


class TestSynctex(TemporaryShowyourworkRepository):
    """Test the paths included in synctex."""

    def check_build(self):
        with gzip.open(self.cwd / "ms.synctex.gz") as f:
            data = f.read().decode("utf-8")

        for line in data.splitlines():
            if not line.startswith("Input:"):
                continue
            path = line.split(":", 2)[-1].strip()
            if not len(path):
                continue
            path = Path(path).relative_to(self.cwd)
            if path.name == "ms.tex":
                assert path.as_posix() == "src/tex/ms.tex"
                return
        raise AssertionError("ms.tex not found in synctex")
