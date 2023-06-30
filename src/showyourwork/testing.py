import filecmp
import hashlib
import os
import platform
import shutil
import subprocess
import sys
import weakref
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory as _TemporaryDirectory
from typing import Any, Generator, Iterable, List, Optional, Union

from showyourwork import cli
from showyourwork.paths import PathLike, find_project_root

# We put all the conda environments in a single directory so that we can reuse
# them between tests
conda_prefix = str(Path().resolve() / ".test" / "conda")


@contextmanager
def cwd(path: PathLike) -> Generator[None, None, None]:
    old_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)


class TemporaryDirectory:
    def __init__(
        self, path: PathLike, args: Iterable[str] = (), force_explicit: bool = False
    ):
        # Snakemake fails when creating conda environments on Windows when executed
        # within a TemporaryDirectory so we wrap it here with a custom
        # implementation on Windows.
        if force_explicit or platform.system() == "Windows":
            m = hashlib.sha1()
            m.update(str(path).encode())
            for arg in args:
                m.update(arg.encode())

            tempdir = Path().resolve() / ".test" / "tmp" / m.hexdigest()
            tempdir.mkdir(parents=True, exist_ok=True)
            self.tempdir = None
            self.name = tempdir

        else:
            self.tempdir = _TemporaryDirectory()
            self.name = Path(self.tempdir.name)

    def cleanup(self) -> None:
        if self.tempdir is None:
            shutil.rmtree(self.name)
        else:
            self.tempdir.cleanup()


class run:
    def __init__(
        self,
        path: PathLike,
        *args: str,
        check_exists: bool = True,
        check_contents: bool = True,
        show_diff: bool = False,
        diff_command: Union[str, Iterable[str]] = ("diff", "-u"),
        expected_dirname: PathLike = "expected",
        **kwargs: Any,
    ):
        self._directory = TemporaryDirectory(path, args)
        self._finalizer = weakref.finalize(self, self._cleanup, self._directory)

        tmpdir = self._directory.name
        test_project_root = Path(path).resolve()

        # Helper function for ignoring the "expected" directory
        def ignore_expected(_: str, names: List[str]) -> List[str]:
            return [
                name for name in names if Path(name).parts[0].startswith("expected")
            ]

        # Copy the test project over to a temporary directory
        shutil.copytree(
            test_project_root, tmpdir, ignore=ignore_expected, dirs_exist_ok=True
        )
        print(list(tmpdir.glob("*")))
        with cwd(tmpdir):
            self.execute(
                *args,
                cwd=tmpdir,
                **kwargs,
            )

        diff_command = [diff_command] if isinstance(diff_command, str) else diff_command
        expected_dir = test_project_root / expected_dirname
        if (check_exists or check_contents) and expected_dir.is_dir():
            diffs: List[Path] = []
            for expected in expected_dir.glob("**/*"):
                # We don't check directories, only files. We can revisit this if
                # necessary.
                if expected.is_dir():
                    continue

                # Construct the path to the expected file in the temporary directory
                subpath = expected.relative_to(test_project_root / expected_dirname)
                observed = tmpdir / subpath

                # Files with the suffix ".exists" are just used to check that
                # the file gets created
                if expected.suffix == ".exists":
                    if not observed.with_suffix("").is_file():
                        raise ValueError(f"{subpath.with_suffix('')} doesn't exist")
                    continue

                if not observed.is_file():
                    if check_exists:
                        raise ValueError(f"{subpath} doesn't exist or is not a file")
                    else:
                        # If we don't care about existence, then we can skip
                        # missing files
                        continue

                if not check_contents:
                    continue

                files_match = filecmp.cmp(expected, observed, shallow=False)
                if not files_match:
                    diffs.append(subpath)

                    if show_diff:
                        diff = _diff(diff_command, expected, observed)
                        print(diff, file=sys.stderr)

            if diffs:
                raise ValueError(
                    "The following files differ from expected contents:\n"
                    + "\n".join(map("- {}".format, diffs))
                )

    def execute(
        self,
        executable: str,
        *args: str,
        cwd: Optional[PathLike] = None,
        **kwargs: Any,
    ) -> None:
        result = subprocess.run(
            [executable, *args],
            check=False,
            capture_output=True,
            text=True,
            cwd=cwd,
            **kwargs,
        )
        if result.returncode:
            raise RuntimeError(
                f"{executable} failed with the following output:\n"
                f"stdout: ===\n{result.stdout}\n===\n\n"
                f"stderr:===\n{result.stderr}\n===\n\n"
            )
        else:
            sys.stderr.write(result.stderr)

    def __enter__(self) -> Path:
        return self._directory.name

    def __exit__(self, *_: Any) -> None:
        self.cleanup()

    @classmethod
    def _cleanup(cls, directory: TemporaryDirectory) -> None:
        directory.cleanup()

    def cleanup(self) -> None:
        if self._finalizer.detach():
            self._cleanup(self._directory)


class run_snakemake(run):
    def __init__(
        self,
        path: PathLike,
        *args: str,
        executable: str = "snakemake",
        check_exists: bool = True,
        check_contents: bool = True,
        show_diff: bool = False,
        diff_command: Union[str, Iterable[str]] = ("diff", "-u"),
        expected_dirname: PathLike = "expected",
        conda_frontend: str = "mamba",
        **kwargs: Any,
    ):
        args = tuple(
            [
                executable,
                "--cores",
                "1",
                "--use-conda",
                "--conda-frontend",
                conda_frontend,
                "--conda-prefix",
                conda_prefix,
            ]
            + list(args)
        )
        super().__init__(
            path,
            *args,
            check_exists=check_exists,
            check_contents=check_contents,
            show_diff=show_diff,
            diff_command=diff_command,
            expected_dirname=expected_dirname,
            **kwargs,
        )


class run_showyourwork(run):
    def __init__(
        self,
        path: PathLike,
        *args: str,
        check_exists: bool = True,
        check_contents: bool = True,
        show_diff: bool = False,
        diff_command: Union[str, Iterable[str]] = ("diff", "-u"),
        expected_dirname: PathLike = "expected",
        configfile: Optional[PathLike] = None,
        cores: str = "1",
        conda_frontend: Optional[str] = "mamba",
    ):
        super().__init__(
            path,
            *args,
            check_exists=check_exists,
            check_contents=check_contents,
            show_diff=show_diff,
            diff_command=diff_command,
            expected_dirname=expected_dirname,
            configfile=configfile,
            cores=cores,
            conda_frontend=conda_frontend,
        )

    def execute(self, *args: Any, **kwargs: Any) -> None:
        kwargs.pop("cwd", None)
        find_project_root.cache_clear()
        cli._build(
            verbose=False,
            snakemake_args=args,
            **kwargs,
        )


def _diff(diff_command: Iterable[str], expected: Path, observed: Path) -> str:
    tempdir = TemporaryDirectory(expected, [str(observed)])
    try:
        a = tempdir.name / "actual" / observed.name
        b = tempdir.name / "expected" / expected.name
        a.parent.mkdir(parents=True, exist_ok=True)
        b.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(observed, a)
        shutil.copy(expected, b)
        result = subprocess.run(
            list(diff_command)
            + [a.relative_to(tempdir.name), b.relative_to(tempdir.name)],
            check=False,
            capture_output=True,
            cwd=tempdir.name,
        )
    finally:
        tempdir.cleanup()
    return result.stdout.decode()
