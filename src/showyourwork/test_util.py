import os
import subprocess
from contextlib import contextmanager
from tempfile import TemporaryDirectory
from typing import Generator, Optional

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
