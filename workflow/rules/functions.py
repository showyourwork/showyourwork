"""
Functions used as inputs / params in various rules.

"""
import json
from pathlib import Path
import warnings
from sphinx_mock import *


def input_class_file(wildcards):
    """
    Return the path to the LaTeX document class file.

    """
    checkpoints.class_name.get(**wildcards)
    with open(relpaths.temp / "class_name", "r") as f:
        folder = f.read().replace("\n", "")
    return posix(abspaths.workflow / "resources" / "classes" / folder / wildcards.file)


def class_files(wildcards):
    """
    Return a list of class files needed to compile the PDF.

    """
    checkpoints.class_name.get(**wildcards)
    with open(relpaths.temp / "class_name", "r") as f:
        folder = f.read().replace("\n", "")
    return [posix(relpaths.tex / file) for file in files.cls.get(folder, [])]


def script_name(wildcards, input):
    """
    Returns the name of the figure script for the `figure` rule.

    Assumes the first python file in the `input` is the one that
    generates the script.

    """
    py_scripts = [file for file in input if file.endswith(".py")]
    return str(Path(py_scripts[0]).relative_to(Path("src") / "figures"))


def figure_script(wildcards):
    """
    Return the figure script that produces `wildcards.figure`.

    """
    checkpoints.script_info.get(**wildcards)
    figure = wildcards.figure
    with open(relpaths.temp / "scripts.json", "r") as f:
        scripts = json.load(f)
    for entry in scripts["figures"].values():
        if figure in entry["files"]:
            script = entry["script"]
            if script != files.unknown:
                if Path(script).exists():
                    return script
                else:
                    raise ShowyourworkException(
                        f"Input script `{script}` for figure `{figure}` does not exist.",
                        rule_name="figure",
                        brief=f"Input script `{script}` for figure `{figure}` does not exist.",
                        context="showyourwork infers the script that generates a given figure "
                        "in your TeX file by inspecting the label inside the figure environment. "
                        r"The convention is that if a figure is labeled `\label{fig:foo}`, "
                        "showyourwork expects there to be a script called `foo.py` in the `src/figures` "
                        "directory that generates the figure(s) in the current figure environment.",
                        delayed=False,
                    )
            else:
                raise ShowyourworkException(
                    f"Unable to generate figure `{figure}`."
                    "Did you forget to provide a custom rule for it?",
                    rule_name="figure",
                    brief=f"Unable to generate figure `{figure}`.",
                    context="showyourwork doesn't know how to generate the "
                    f"figure `{figure}` in your TeX file. "
                    "Usually, showyourwork is able to infer the script that "
                    r"generates a given figure by parsing \includegraphics "
                    r"and \label commands in your TeX file, provided they are "
                    "inside a figure environment. But if a figure is included "
                    "outside of a figure environment, showyourwork can't "
                    "automatically infer its parent script. In such cases, you "
                    "must provide a custom rule in your Snakefile to generate "
                    "the figure; see the documentation for details.",
                    delayed=False,
                )
    raise ShowyourworkException(
        f"Input script not found for figure `{figure}`.",
        rule_name="figure",
        brief=f"Input script not found for figure `{figure}`.",
        context="showyourwork infers the script that generates a given figure "
        "in your TeX file by inspecting the label inside the figure environment. "
        r"The convention is that if a figure is labeled `\label{fig:foo}`, "
        "showyourwork expects there to be a script called `foo.py` in the `src/figures` "
        "directory that generates the figure(s) in the current figure environment.",
        delayed=False,
    )


def figure_script_dependencies(wildcards):
    """
    Return user-specified dependencies of the current figure script.

    """
    script = Path(figure_script(wildcards)).name
    deps = []
    for dep in config["dependencies"].get(script, []):
        if type(dep) is OrderedDict:
            dep = list(dep)[0]
        deps.append(str(relpaths.src / dep))
    return deps


def script_dependencies(wildcards):
    """
    Return user-specified dependencies of the current script.

    """
    script = zenodo.script[wildcards.dependency]
    deps = [str(relpaths.src / script)]
    for dep in config["dependencies"].get(script, []):
        if type(dep) is OrderedDict:
            dep = list(dep)[0]
        deps.append(str(relpaths.src / dep))
    return deps


def figures(wildcards):
    """
    Return all the figure files required by the manuscript.

    """
    checkpoints.script_info.get(**wildcards)
    figures = []
    with open(relpaths.temp / "scripts.json", "r") as f:
        scripts = json.load(f)
    for entry in scripts["figures"].values():
        figures += entry["files"]
    return figures


def check_figure_format(figure):
    """
    Check that all figures are declared correctly in `tex/ms.tex`
    so we can parse them corresponding XML tree.

    """
    # Get all figure elements
    elements = list(figure)
    captions = figure.findall("CAPTION")
    labels = figure.findall("LABEL")

    # Check that figure labels aren't nested inside captions
    for caption in captions:
        caption_labels = caption.findall("LABEL")
        if len(caption_labels):
            raise ValueError(
                "Label `{}` should not be nested within the figure caption".format(
                    caption_labels[0].text
                )
            )

    # The label must always come after the figure caption
    if len(captions):

        # Index of last caption
        for caption_idx, element in enumerate(elements[::-1]):
            if element.tag == "CAPTION":
                break
        caption_idx = len(elements) - 1 - caption_idx

        if len(labels):

            # Index of first label
            for label_idx, element in enumerate(elements):
                if element.tag == "LABEL":
                    break

            if label_idx < caption_idx:
                raise ValueError(
                    "Figure label `{}` must come after the caption.".format(
                        (labels)[0].text
                    )
                )

    # Check that there is exactly one label
    if len(labels) >= 2:
        raise ValueError(
            "A figure has multiple labels: `{}`".format(
                ", ".join(label.text for label in labels)
            )
        )
    elif len(labels) == 0:
        # Unable to infer script
        warnings.warn("There is a figure without a label.")