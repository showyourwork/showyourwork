import subprocess
from abc import ABC, abstractmethod
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Generator, Iterable, Optional, Union

from showyourwork.paths import PathLike


class vcs(ABC):
    @abstractmethod
    def root(self) -> Path:
        ...

    @abstractmethod
    def branch(self) -> str:
        ...

    @abstractmethod
    def hash(self) -> str:
        ...

    @abstractmethod
    def tag(self) -> Optional[str]:
        ...


class git(vcs):
    def __init__(self, executable: Optional[PathLike] = None):
        self.executable = "git" if executable is None else str(executable)

    def git(
        self, args: Union[Iterable[str], str], **kwargs: Any
    ) -> subprocess.CompletedProcess[str]:
        kwargs["check"] = kwargs.get("check", True)
        kwargs["text"] = kwargs.get("text", True)
        kwargs["capture_input"] = kwargs.get("capture_input", True)
        if isinstance(args, str):
            args = f"{self.executable} {args}"
            kwargs["shell"] = kwargs.get("shell", True)
        else:
            args = [self.executable, *args]
            kwargs["shell"] = kwargs.get("shell", False)
        return subprocess.run(args, **kwargs)

    def root(self) -> Path:
        return Path(self.git(["rev-parse", "--show-toplevel"]).stdout.strip())

    def branch(self) -> str:
        return self.git(["rev-parse", "--abbrev-ref", "HEAD"]).stdout.strip()

    def hash(self) -> str:
        return self.git(["rev-parse", "HEAD"]).stdout.strip()

    def tag(self) -> Optional[str]:
        r = self.git(["describe", "--exact-match", "--tags", "HEAD"], check=False)
        if r.returncode == 0:
            return r.stdout.strip()
        return None

    @contextmanager
    def checkout_temp(self, branch: Optional[str] = None) -> Generator[str, None, None]:
        args = ["clone", "--depth", "1"]
        if branch is not None:
            args += ["--branch", branch]
        args += [str(self.root())]
        with TemporaryDirectory() as d:
            self.git(args + [d])
            yield d
