"""
Helper functions for Snakemake.

Note that this file is imported directly into the Snakefile,
so all paths are relative to the root directory.

"""
import shutil
from pathlib import Path
import sys

sys.path.insert(1, ".cortex")
from cortex import (
    get_scripts,
    gen_pdf,
    get_user_metadata,
    get_script_metadata,
    get_metadata,
)


def scripts(wildcards, tags=["figures", "tests"]):
    checkpoints.check_scripts.get(**wildcards)
    return list(get_scripts(tags=tags).keys())


def figures(wildcards, flat=True):
    checkpoints.check_scripts.get(**wildcards)
    figures = list(get_scripts(tags=["figures"]).values())
    if flat:
        return [lst for entry in figures for lst in entry]
    else:
        return figures


def figure_script(wildcards):
    checkpoints.check_scripts.get(**wildcards)
    figure = wildcards.figure
    scripts = get_scripts(tags=["figures"])
    for script, files in scripts.items():
        if figure in files:
            return script
    raise ValueError(
        "Input script not found for output figure `{}`.".format(figure)
    )


def figure_script_base_name(wildcards):
    return Path(figure_script(wildcards)).name


def figure_cache(wildcards, output=None):
    return Path(".cortex") / "data" / Path(output[0]).name


def figure_other(wildcards, output=None):
    checkpoints.check_scripts.get(**wildcards)
    all_output = get_scripts(tags=["figures"])[figure_script(wildcards)]
    return [Path(file) for file in all_output if file != output[0]]


def run_figure(params):
    if params.cache.exists():
        shutil.move(params.cache, "figures")
    else:
        shell("cd figures && python {params.script}")
        for file in params.other:
            if file.exists():
                shutil.move(file, Path(".cortex") / "data")
