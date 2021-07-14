"""
Helper functions for Snakemake.

Note that this file is imported directly into the Snakefile,
so all paths are relative to the root directory.

"""
import shutil
from pathlib import Path
import sys


sys.path.insert(1, ".cortex")
import cortex

cortex.config = config

from cortex import (
    get_scripts,
    gen_pdf,
    get_user_metadata,
    get_script_metadata,
    get_metadata,
)


def scripts(wildcards, tags=["figures", "tests"]):
    checkpoints.script_meta.get(**wildcards)
    return list(get_scripts(tags=tags).keys())


def figures(wildcards):
    checkpoints.script_meta.get(**wildcards)
    figures = list(get_scripts(tags=["figures"]).values())
    return [lst for entry in figures for lst in entry]


def test_results(wildcards):
    checkpoints.script_meta.get(**wildcards)
    tests = list(get_scripts(tags=["tests"]).values())
    return [lst for entry in tests for lst in entry]


def figure_script(wildcards):
    checkpoints.script_meta.get(**wildcards)
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
    checkpoints.script_meta.get(**wildcards)
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


figure_wildcards = "figures/(.*?)\.{}".format(
    "|".join(
        config.get(
            "figure_extensions",
            ["pdf", "png", "epsfont", "jpg", "jpeg", "gif"],
        )
    )
)


def run_test(input, output):
    passed = r"\\\def\\\ctxCurrentTestBadge{\\\color{ctxTestPassed}\\\faCheck}"
    failed = r"\\\def\\\ctxCurrentTestBadge{\\\color{ctxTestFailed}\\\faTimes}"
    shell(
        "{{ py.test {input[0]} &>/dev/null && echo {passed} || echo {failed} ; }} > {output[0]}"
    )


test_wildcards = "tests/test_(.*?)\.py"
