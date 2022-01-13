import subprocess
import shutil
import sys
import json
import re
import warnings
from pathlib import Path
from xml.etree.ElementTree import parse as ParseXMLTree

# Import utils
sys.path.insert(1, snakemake.config["workflow_path"])
from utils import paths


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
                raise ValueError(
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
                raise ValueError(
                    "Command \marginicon must always come before the figure label."
                )

    # Check that there is exactly one label
    if len(labels) >= 2:
        raise ValueError(
            "A figure has multiple labels: `{}`".format(
                ", ".join(label.text for label in labels)
            )
        )
    elif len(labels) == 0:
        if len(figure.findall("LABELSTAR")) == 0:
            raise ValueError("There is a figure without a label.")


def get_xml_tree():
    """"""
    # Parameters
    tex_files_in = snakemake.config["tex_files_in"]
    ms_name = snakemake.config["ms_name"]
    ms = paths.tex / f"{ms_name}.tex"
    xmlfile = paths.temp / "showyourwork.xml"

    # Copy over TeX auxiliaries
    for file in tex_files_in:
        src = Path(file)
        dst = paths.tex / src.name
        if not dst.exists():
            shutil.copy(str(src), str(dst))

    # Run tectonic
    result = subprocess.run(
        ["tectonic", "--chatter", "minimal", "-r", "0", "-o", str(paths.temp), ms],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        # TODO: Log these
        # TODO: Keep all logs
        print(result.stdout.decode("utf-8"))
        print(result.stderr.decode("utf-8"))
        raise ValueError

    # Add <HTML></HTML> tags to the XML file
    try:
        with open(xmlfile, "r") as f:
            contents = f.read()
    except FileNotFoundError:
        contents = ""
    contents = "<HTML>\n" + contents + "</HTML>"
    with open(xmlfile, "w") as f:
        print(contents, file=f)

    # Load the XML tree
    return ParseXMLTree(paths.temp / "showyourwork.xml").getroot()


def get_json_tree(xml_tree):
    """"""
    # Parse the \graphicspath command
    # Note that if there are multiple graphicspath calls, only the first one
    # is read. Same for multiple directories within a graphicspath call.
    try:
        graphicspath = xml_tree.findall("GRAPHICSPATH")
        if len(graphicspath) == 0:
            graphicspath = Path(".")
        else:
            graphicspath = re.findall("\{(.*?)\}", graphicspath[-1].text)[0]
            graphicspath = Path(graphicspath)
    except:
        warnings.warn("Unable to parse the `graphicspath` command. Ignoring...")
        graphicspath = Path(".")

    # Parse labeled graphics inside `figure` environments
    figures = {}
    for figure in xml_tree.findall("FIGURE"):

        # Ensure the figure environment conforms to the standard
        check_figure_format(figure)

        # Get the figure \label
        labels = figure.findall("LABEL")
        label_stars = figure.findall("LABELSTAR")

        # Infer the figure script from the label
        if len(labels) and labels[0].text is not None:

            # User defined a \label{fig:XXX}
            label = labels[0].text

            # Infer the full path to the script, searching for
            # files with any of the user-defined extensions
            # (and settling on the first match)
            for ext in snakemake.config["script_extensions"]:
                script = paths.figure_scripts / f"{label}.{ext}"
                if script.exists():
                    script = str(script.relative_to(paths.user))
                    break
            else:
                # Fallback to ".py"; we'll catch the error later!
                script = str(
                    (paths.figure_scripts / f"{label}.py").relative_to(paths.user)
                )

        elif len(label_stars) and label_stars[0].text is not None:

            # User defined a \label{fig*:XXX}
            label = labels[0].text

            # There is no script associated with this figure
            script = None

        else:

            raise ValueError("There is a figure without a label.")

        # Find all graphics included in this figure environment
        graphics = [
            str(
                (paths.tex / graphicspath / graphic.text)
                .resolve()
                .relative_to(paths.user)
            )
            for graphic in figure.findall("GRAPHICS")
        ]

        # TODO: Collect associated datasets
        datasets = []

        figures[label] = {"script": script, "graphics": graphics, "datasets": datasets}

    # Parse free-floating graphics
    graphics = [
        str((paths.tex / graphicspath / graphic.text).resolve().relative_to(paths.user))
        for graphic in xml_tree.findall("GRAPHICS")
    ]
    figures["free-floating"] = {"script": None, "graphics": graphics, "datasets": []}

    # The full tree (someday we'll have equations in here, too)
    tree = {"figures": figures}

    return tree


# Get the article tree
xml_tree = get_xml_tree()
snakemake.config["tree"] = get_json_tree(xml_tree)


# Save the config file
with open(paths.temp / "config.json", "w") as f:
    print(json.dumps(snakemake.config, indent=4), file=f)