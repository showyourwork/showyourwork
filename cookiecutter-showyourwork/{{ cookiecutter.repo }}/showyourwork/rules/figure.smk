import json
from pathlib import Path


def script_name(wildcards, input):
    """
    Returns the name of the figure script for the `figure` rule.

    Assumes the first python file in the `input` is the one that
    generates the script.

    """
    py_scripts = [file for file in input if file.endswith(".py")]
    return Path(py_scripts[0]).name


def figure_script(wildcards):
    """
    Return the figure script that produces `wildcards.figure`.

    """
    checkpoints.script_info.get(**wildcards)
    figure = wildcards.figure
    with open(TEMP / "scripts.json", "r") as f:
        scripts = json.load(f)
    for entry in scripts["figures"].values():
        if figure in entry["files"]:
            return entry["script"]
    raise ValueError(
        "Input script not found for output figure `{}`.".format(figure)
    )


def figures(wildcards):
    """
    Return all the figure files required by the manuscript.

    """
    checkpoints.script_info.get(**wildcards)
    figures = []
    with open(TEMP / "scripts.json", "r") as f:
        scripts = json.load(f)
    for entry in scripts["figures"].values():
        figures += entry["files"]
    return figures


rule figure:
    message:
        "Generating figure `{output}`..."
    input:
        figure_script,
        "environment.yml",
    output:
        "{figure}",
    wildcard_constraints:
        figure="figures/(.*?)\.{}".format("|".join(figexts)),
    params:
        script_name=script_name,
    conda:
        POSIX(USER / "environment.yml")
    shell:
        "cd figures && python {params.script_name}"
