import json
import os
from xml.etree.ElementTree import parse as ParseXMLTree


localrules: script_info


checkpoint script_info:
    """
    Build the figure dependency tree.

    """
    message:
        "Building figure dependency tree..."
    input:
        posix(relpaths.temp / "showyourwork.xml"),
        files.dot_zenodo
    output:
        posix(relpaths.temp / "scripts.json"),
    run:
        # Load the XML tree
        root = ParseXMLTree(relpaths.temp / "showyourwork.xml").getroot()

        # Parse graphics inside `figure` environments
        figures = {}
        for figure in root.findall("FIGURE"):
            check_figure_format(figure)
            labels = figure.findall("LABEL")
            if len(labels):
                label = labels[0].text
                if label is not None:
                    script = posix(relpaths.figures / "{}.py".format(label))
                    files = []
                    for graphic in figure.findall("GRAPHICS"):
                        if graphic.text.startswith("figures/"):
                            files.append(f"src/{graphic.text}")
                        elif graphic.text.startswith("static/"):
                            continue
                        elif graphic.text.lower() in files.special_figures:
                            continue
                        else:
                            warnings.warn(
                                f"Figure `{graphic.text}` must be in either the `src/figures` or `src/static` folders."
                            )
                    if len(files):
                        datasets = [
                            dataset for dataset in 
                            config["figure_dependencies"].get(
                                "{}.py".format(label), []
                            ) if type(dataset) is str and 
                            dataset.endswith(".zenodo")
                        ]
                        if label in figures:
                            figures[label]["files"] += files
                            figures[label]["datasets"] += datasets
                        else:
                            figures[label] = {
                                "script": script, 
                                "files": files, 
                                "datasets": datasets
                            }

        # Parse free-floating graphics
        files = []
        for graphic in root.findall("GRAPHICS"):
            if graphic.text.startswith("figures/"):
                files.append(f"src/{graphic.text}")
            elif graphic.text.startswith("static/"):
                continue
            elif graphic.text.lower() in files.special_figures:
                continue
            else:
                warnings.warn(
                    f"Figure `{graphic.text}` must be in either the `src/figures` or `src/static` folders."
                )
        if len(files):
            figures["unknown"] = {
                "script": files.unknown,
                "files": files,
                "datasets": []
            }

        # Store as JSON
        scripts = {"figures": figures}
        with open(relpaths.temp / "scripts.json", "w") as f:
            print(json.dumps(scripts, indent=4), file=f)
