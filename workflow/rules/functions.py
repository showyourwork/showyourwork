"""
Functions used as inputs / params in various rules.

"""
import json
from pathlib import Path
import warnings
from sphinx_mock import *


class File:
    """
    Object used to resolve placeholder strings in shell commands.

    """

    def __init__(self, file, **kwargs):
        self.file = file
        self.path = str(Path(file).parents[0])
        self.name = Path(file).name
        self.extension = Path(file).suffix.split(".")[-1]
        self.__dict__.update(**kwargs)

    def __repr__(self):
        return str(self.file)


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
    Returns the name of the script for the `figure` rule.

    Assumes the first *script* in the `input` is the one that
    generates the output.

    """
    # Get all files with these extensions
    scripts = []
    for file in input:
        for ext in files.script_extensions:
            if file.endswith(ext):
                scripts.append(file)
    return str(Path(scripts[0]).relative_to(Path("src") / "figures"))


def shell_cmd(wildcards, input, output):
    """
    Returns the shell command to produce an output from a script.

    """
    # Get the full path to the output file
    output_file = output[0]

    # Get the full path to the input script
    if hasattr(wildcards, "figure"):
        # If this is being run from the `figure` rule, the
        # name of the script can be automatically inferred
        input_script = figure_script(wildcards)
    elif hasattr(wildcards, "dependency"):
        # If this is being run from the `upload` rule, we
        # can also easily infer the name of the script
        input_script = zenodo.script[wildcards.dependency]
    else:
        # Assume the script that generates the output is the
        # first script in the list of inputs
        try:
            for file in input:
                for ext in files.script_extensions:
                    if file.endswith(ext):
                        input_script = file
                        raise StopIteration
        except StopIteration:
            pass
        else:
            raise ShowyourworkException(
                f"Cannot determine input script for output file {output_file}.",
                brief=f"Unable to determine input script.",
                delayed=False,
            )

    # If the user specifies a custom Snakemake rule, there may not
    # be a script and that's ok
    if input_script == files.unknown:
        return ""

    # Get the command to run this kind of script
    ext = Path(input_script).suffix.split(".")[-1]
    cmd = config["scripts"].get(ext, None)

    if cmd is None:
        raise ShowyourworkException(
            f"Unknown script extension: `{ext}`.",
            brief=f"Unknown script extension: `{ext}`.",
            context="showyourwork does not know how to generate output from "
            f"scripts ending in `{ext}`. Please provide instructions "
            "in the `showyourwork.yml` config file; see the docs for details.",
            delayed=False,
        )

    # Fill in the placeholders
    return cmd.format(script=File(input_script), output=File(output_file))


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
                        "directory that generates the figure(s) in the current figure environment. "
                        "Note that showyourwork supports figure scripts with other extensions, but in "
                        "those cases users must explicitly tell showyourwork how to generate figures "
                        "from those scripts by providing a shell command in `showyourwork.yml`. Please "
                        "see the docs for details. If you're still getting this error message after "
                        "fixing the issue, try deleting the `.showyourwork` folder to clear the build "
                        "cache.",
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
                    delayed=True,
                )
    raise ShowyourworkException(
        f"Input script not found for figure `{figure}`.",
        rule_name="figure",
        brief=f"Input script not found for figure `{figure}`.",
        context="showyourwork infers the script that generates a given figure "
        "in your TeX file by inspecting the label inside the figure environment. "
        r"The convention is that if a figure is labeled `\label{fig:foo}`, "
        "showyourwork expects there to be a script called `foo.py` in the `src/figures` "
        "directory that generates the figure(s) in the current figure environment. "
        "Note that showyourwork supports figure scripts with other extensions, but in "
        "those cases users must explicitly tell showyourwork how to generate figures "
        "from those scripts by providing a shell command in `showyourwork.yml`. Please "
        "see the docs for details. If you're still getting this error message after "
        "fixing the issue, try deleting the `.showyourwork` folder to clear the build "
        "cache.",
        delayed=False,
    )


def script_dependencies(wildcards):
    """
    Return user-specified dependencies of a figure or a Zenodo dataset.

    """
    # Get the full path to the input script
    if hasattr(wildcards, "figure"):
        # This is a figure dependency
        script = figure_script(wildcards)
        deps = []
    else:
        # This is a Zenodo dependency
        script = zenodo.script[wildcards.dependency]
        deps = [str(script)]  # this is the script that generates the dataset
    for dep in config["dependencies"].get(script, []):
        if type(dep) is OrderedDict:
            dep = list(dep)[0]
        deps.append(str(dep))
    return deps


def dot_zenodo_files(wildcards):
    checkpoints.script_info.get(**wildcards)
    datasets = []
    with open(relpaths.temp / "scripts.json", "r") as f:
        scripts = json.load(f)
    for _, entry in scripts["figures"].items():
        datasets += entry["datasets"]
    return datasets


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
    # Any marginicons must always come before the label
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
                raise ShowyourworkException(
                    "Figure label `{}` must come after the caption.".format(
                        (labels)[0].text
                    )
                )

            # Index of last marginicon
            for marginicon_idx, element in enumerate(elements):
                if element.tag == "MARGINICON":
                    break
            else:
                marginicon_idx = 0

            if marginicon_idx > label_idx:
                raise ShowyourworkException(
                    "Command \marginicon must come before the figure label."
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
        # Only warn the user if they didn't provide a \label{fig*:...}
        if len(figure.findall("LABELSTAR")) == 0:
            warnings.warn("There is a figure without a label.")