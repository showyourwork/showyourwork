import os
import subprocess
import sys

from ... import paths


def run_snakemake(
    snakefile,
    run_type=None,
    cores=1,
    conda_frontend="conda",
    extra_args=[],
    check=True,
):
    env = dict(os.environ)
    env["SNAKEMAKE_OUTPUT_CACHE"] = paths.user().cache.as_posix()
    if run_type is not None:
        env["SNAKEMAKE_RUN_TYPE"] = run_type
    cmd = [
        "snakemake",
        "--cores",
        f"{cores}",
        "--use-conda",
        "--conda-frontend",
        conda_frontend,
        "--reason",
        "--cache",
        "-s",
        snakefile,
    ] + list(extra_args)
    result = subprocess.run(cmd, env=env, check=False)
    if check and result.returncode:
        sys.exit(result.returncode)
    return result.returncode
