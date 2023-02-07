import os
import subprocess
import sys
from typing import Iterable, Optional

import click

from showyourwork import paths
from showyourwork.version import __version__


@click.group()
@click.version_option(
    __version__,
    "--version",
    "-v",
    package_name="showyourwork",
    message="%(version)s",
)
def main() -> None:
    """Easily build open-source, reproducible scientific articles."""
    pass


@main.command(
    context_settings=dict(
        ignore_unknown_options=True,
    )
)
@click.option(
    "--configfile",
    type=click.Path(exists=True),
    help="A showyourwork configuration file",
)
@click.option(
    "-c",
    "--cores",
    default="1",
    help="Number of cores to use; passed to snakemake",
)
@click.option(
    "--conda-frontend",
    default="conda",
    help="The conda frontend to use; passed to snakemake",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def build(
    configfile: Optional[paths.PathLike],
    cores: str,
    conda_frontend: str,
    snakemake_args: Iterable[str],
) -> None:
    """Build an article in the current working directory."""
    run_snakemake(
        paths.package_data("showyourwork", "workflow", "Snakefile"),
        config_file=configfile,
        cores=cores,
        conda_frontend=conda_frontend,
        check=True,
        extra_args=snakemake_args,
    )


def run_snakemake(
    snakefile: paths.PathLike,
    config_file: Optional[paths.PathLike] = None,
    cores: str = "1",
    conda_frontend: str = "conda",
    check: bool = True,
    extra_args: Iterable[str] = (),
) -> int:
    cwd = paths.find_project_root()
    if config_file is None:
        config_file = cwd / "showyourwork.yml"
        if not config_file.is_file():
            config_file = cwd / "showyourwork.yaml"
        if not config_file.is_file():
            raise RuntimeError(
                "No config file found in project root. Please specify a configuration "
                "file using the '--configfile' command line argument."
            )

    env = dict(os.environ)
    cmd = [
        "snakemake",
        "--cores",
        f"{cores}",
        "--use-conda",
        "--conda-frontend",
        conda_frontend,
        "--reason",
        "--configfile",
        config_file,
        "-s",
        snakefile,
    ] + list(extra_args)
    result = subprocess.run(cmd, env=env, check=False, cwd=cwd)
    if check and result.returncode:
        sys.exit(result.returncode)
    return result.returncode
