import difflib
import os
import shutil
import subprocess
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Generator, Iterable, List, Optional

from showyourwork import cli
from showyourwork.paths import PathLike, find_project_root

conda_temporary_directory = TemporaryDirectory(prefix="showyourwork-conda-")


@contextmanager
def temporary_project(config: str = "config-version: 2") -> Generator[str, None, None]:
    old_cwd = os.getcwd()
    try:
        find_project_root.cache_clear()
        with TemporaryDirectory() as d:
            open(f"{d}/showyourwork.yml", "w").write(config)
            os.chdir(d)
            yield d
    finally:
        os.chdir(old_cwd)


def run_snakemake(
    snakefile: PathLike,
    targets: list[str],
    conda_frontend: str = "conda",
    cwd: Optional[PathLike] = None,
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        [
            "snakemake",
            "--cores",
            "1",
            "--use-conda",
            "--conda-frontend",
            conda_frontend,
            "--conda-prefix",
            conda_temporary_directory.name,
            "--snakefile",
            str(snakefile),
            *targets,
        ],
        check=False,
        capture_output=True,
        text=True,
        cwd=cwd,
    )
    if result.returncode:
        raise RuntimeError(
            "Snakemake failed with the following output:\n"
            f"stdout: ===\n{result.stdout}\n===\n\n"
            f"stderr:===\n{result.stderr}\n===\n\n"
        )
    return result


@contextmanager
def cwd(path: PathLike) -> Generator[None, None, None]:
    old_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)


def run(
    path: PathLike,
    check_exists: bool = True,
    check_contents: bool = True,
    configfile: Optional[PathLike] = None,
    cores: str = "1",
    conda_frontend: Optional[str] = "mamba",
    snakemake_args: Iterable[str] = (),
) -> None:
    # We need to copy the full project even if the path we've been handed is a
    # subdirectory
    find_project_root.cache_clear()
    test_project_root = find_project_root(path)
    find_project_root.cache_clear()

    # If the test project doesn't have a showyourwork.yml file, then this will
    # sometimes find the top repo directory. We don't want that.
    assert not (test_project_root / ".git").is_dir()

    # Helper function for ignoring the "expected" directory
    def ignore_expected(_: str, names: List[str]) -> List[str]:
        return [name for name in names if Path(name).parts[0] == "expected"]

    # Copy the test project over to a temporary directory
    with TemporaryDirectory() as tmpdir:
        shutil.copytree(
            test_project_root, tmpdir, ignore=ignore_expected, dirs_exist_ok=True
        )
        with cwd(tmpdir):
            cli._build(
                configfile=configfile,
                cores=cores,
                conda_frontend=conda_frontend,
                snakemake_args=snakemake_args,
            )

        if check_exists or check_contents:
            for expected in (test_project_root / "expected").glob("**/*"):
                # We don't check directories, only files. We can revisit this if
                # necessary.
                if expected.is_dir():
                    continue

                # Construct the path to the expected file in the temporary directory
                subpath = expected.relative_to(test_project_root / "expected")
                observed = Path(tmpdir) / subpath

                if not observed.is_file():
                    if check_exists:
                        raise ValueError(f"{subpath} doesn't exist or is not a file")
                    else:
                        # If we don't care about existence, then we can skip
                        # missing files
                        continue

                if not check_contents:
                    continue

                # Compare the contents of the expected and observed files.
                # TODO(dfm): We're assuming that these are plaintext files, but
                # do we want to handle binary files too?
                expected_contents = expected.read_text().strip()
                observed_contents = observed.read_text().strip()
                if expected_contents != observed_contents:
                    diff = difflib.unified_diff(
                        observed_contents.splitlines(keepends=True),
                        expected_contents.splitlines(keepends=True),
                        tofile=f"{expected.relative_to(test_project_root)}",
                        fromfile=f"actual/{subpath}",
                    )
                    raise ValueError(
                        f"{subpath} doesn't match expected output; diff:\n\n"
                        + "".join(diff)
                    )
