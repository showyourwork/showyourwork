import glob
from pathlib import Path
import json
import os
import shutil
import sys
sys.path.insert(1, ".cortex")
import cortex


# Config
configfile: "config.yml"
debug = config["debug"]
quiet = config["quiet"]
ext = "|".join(config["figure_extensions"])


def scripts(wildcards, tags=["figures", "tests"]):
    checkpoints.check_scripts.get(**wildcards)
    return list(cortex.get_scripts(tags=tags).keys())


def figures(wildcards, flat=True):
    checkpoints.check_scripts.get(**wildcards)
    figures = list(cortex.get_scripts(tags=["figures"]).values())
    if flat:
        return [lst for entry in figures for lst in entry]
    else:
        return figures


rule pdf:
    input: 
        "tex/ms.tex", 
        "tex/bib.bib",
        ".cortex/data/meta.json", 
        "figures/fractal_mandelbrot.pdf", 
        "figures/fractal_koch.pdf", 
        "figures/heart.pdf"
    output: 
        "ms.pdf"
    run: 
        cortex.gen_pdf()


rule user_meta:
    output: 
        ".cortex/data/user.json"
    run: 
        cortex.get_user_metadata(clobber=True)


checkpoint check_scripts:
    input: 
        "tex/ms.tex"
    output: 
        ".cortex/data/scripts.json"
    run: 
        cortex.get_script_metadata(clobber=True)


def figure_script(wildcards):
    checkpoints.check_scripts.get(**wildcards)
    figure = wildcards.figure
    scripts = cortex.get_scripts(tags=["figures"])
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


def figure_other(wildcards, output=None, input=None):
    # TODO: Checkpointing
    all_output = cortex.get_scripts(tags=["figures"])[input[0]]
    return [Path(file) for file in all_output if file != output[0]]


rule figure:
    input: 
        figure_script
    output: 
        "{figure}"
    wildcard_constraints:
        figure="figures/(.*?)\.{}".format(ext)
    params:
        script=figure_script_base_name,
        cache=figure_cache,
        other=figure_other
    run:
        if params.cache.exists():
            shutil.move(params.cache, "figures")
        else:
            shell("cd figures && python {params.script}")
            for file in params.other:
                if file.exists():
                    shutil.move(file, Path(".cortex") / "data")



rule meta:
    input: 
        ".cortex/data/user.json", scripts
    output: 
        ".cortex/data/meta.json"
    run: 
        cortex.get_metadata(clobber=True)