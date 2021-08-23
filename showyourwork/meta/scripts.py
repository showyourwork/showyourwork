from ..constants import *
from ..utils import make_pdf, save_json
from pathlib import Path
import subprocess
import json
from xml.etree.ElementTree import parse as ParseXMLTree
import numpy as np
import warnings
import os


__all__ = ["get_script_metadata", "get_script_status"]


def check_figure_format(figure):

    # Get all figure elements
    elements = list(figure)
    captions = figure.findall("CAPTION")
    labels = figure.findall("LABEL") + figure.findall("LABEL_")

    # Check that figure labels aren't nested inside captions
    for caption in captions:
        caption_labels = caption.findall("LABEL") + caption.findall("LABEL_")
        if len(caption_labels):
            raise ValueError(
                "Label `{}` should not be nested within the figure caption".format(
                    caption_labels[0].text
                )
            )

    # The label must always come after the figure caption
    if len(captions):

        # Index of last caption
        caption_idx = (
            len(elements)
            - 1
            - np.argmax(
                [element.tag == "CAPTION" for element in elements[::-1]]
            )
        )

        if len(labels):

            # Index of first label
            label_idx = np.argmax(
                [
                    element.tag == "LABEL" or element.tag == "LABEL_"
                    for element in elements
                ]
            )

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


def get_figure_label(figure):
    labels = figure.findall("LABEL")
    if len(labels):
        return labels[0].text
    else:
        return None


def get_figure_script(label):
    return str(Path("figures") / "{}.py".format(label))


def get_figure_files(figure):
    files = []
    for graphic in figure.findall("GRAPHICS"):
        if graphic.text.startswith("figures/"):
            files.append(graphic.text.replace("/", os.sep))
        else:
            files.append(str(Path("figures") / graphic.text))
    return files


def get_script_metadata(clobber=True, verbose=False):

    if clobber or not (TEMP / "scripts.json").exists():

        # Generate the XML file
        make_pdf(
            tmpdir=TEMP / "tree",
            publish=False,
            verbose=verbose,
            tectonic_args=["--keep-logs", "--keep-intermediates", "-r", "0"],
            gen_tree=True,
        )
        root = ParseXMLTree(TEMP / "tree" / "showyourwork.xml").getroot()

        # Parse graphics inside `figure` environments
        figures = {}
        for figure in root.findall("FIGURE"):
            check_figure_format(figure)
            label = get_figure_label(figure)
            if label is not None:
                script = get_figure_script(label)
                files = get_figure_files(figure)
                figures[label] = {"script": script, "files": files}

        # Parse free-floating graphics
        figures["unknown"] = {
            "script": UNKNOWN_SCRIPT,
            "files": get_figure_files(root),
        }

        # Store as JSON
        scripts = {"figures": figures}
        save_json(scripts, TEMP / "scripts.json")

        return scripts

    else:

        with open(TEMP / "scripts.json", "r") as f:
            scripts = json.load(f)

        return scripts


def get_script_status(script):
    """
    Return an error code corresponding to the git status of a given script.

    """
    # Check if the file exists
    script = Path(script)
    if not script.exists():
        error = ScriptDoesNotExist
    else:
        # Check if the file is version controlled
        try:
            status = (
                subprocess.check_output(
                    ["git", "ls-files", "--error-unmatch", script],
                    cwd=USER,
                    stderr=subprocess.DEVNULL,
                )
                .decode()
                .replace("\n", "")
            )
        except Exception as e:
            # File is not version controlled
            error = ScriptNotVersionControlled
        else:
            # Check if the file has uncommitted changes
            try:
                status = (
                    subprocess.check_output(
                        ["git", "status", "-s", script],
                        cwd=USER,
                        stderr=subprocess.DEVNULL,
                    )
                    .decode()
                    .replace("\n", "")
                )
                if len(status):
                    raise Exception("Uncommitted changes!")
            except Exception as e:
                # Uncommited changes
                error = ScriptHasUncommittedChanges
            else:
                # File is good!
                error = ScriptUpToDate
    return error
