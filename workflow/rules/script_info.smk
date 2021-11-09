import json
import os
import re
from xml.etree.ElementTree import parse as ParseXMLTree


localrules: script_info


checkpoint script_info:
    """
    Builds the figure dependency tree from the XML file
    generated in the ``xml`` rule. Stores it in a JSON
    file in the temporary ``showyourwork`` directory.

    """
    message:
        "Building figure dependency tree..."
    input:
        posix(relpaths.temp / "showyourwork.xml"),
        [posix(f) for f in files.dot_zenodo]
    output:
        posix(relpaths.temp / "scripts.json"),
    run:
        # Load the XML tree
        root = ParseXMLTree(relpaths.temp / "showyourwork.xml").getroot()

        # Parse graphicspath command
        # NOTE: If there are multiple graphicspath calls, only the first one
        # is read. Same for multiple directories within a graphicspath call.
        try:
            graphicspath = root.findall("GRAPHICSPATH")
            if len(graphicspath) == 0:
                graphicspath = relpaths.figures
            else:
                graphicspath = re.findall("\{(.*?)\}", graphicspath[-1].text)[0]
                graphicspath = Path(graphicspath)
        except:
            warnings.warn("Unable to parse the `graphicspath` command. Ignoring...")
            graphicspath = Path(".")

        # Parse graphics inside `figure` environments
        figures = {}
        for figure in root.findall("FIGURE"):
            check_figure_format(figure)
            labels = figure.findall("LABEL")
            if len(labels):
                label = labels[0].text
                if label is not None:
                    for ext in files.script_extensions:
                        script = posix(relpaths.figures / "{}.{}".format(label, ext))
                        if Path(script).exists():
                            break
                    else:
                        # Fallback to ".py"; we'll catch the error later!
                        script = posix(relpaths.figures / "{}.py".format(label))
                    filenames = []
                    for graphic in figure.findall("GRAPHICS"):
                        src = relpaths.src.absolute()
                        path = (src / graphicspath / graphic.text).resolve().relative_to(src)
                        if path.parts[0] == "figures":
                            filenames.append(str(relpaths.src / path))
                        elif path.parts[0] == "static":
                            continue
                        elif path.name in files.special_figures:
                            continue
                        else:
                            warnings.warn(
                                f"Figure `{graphic.text}` must be in either the `src/figures` or `src/static` folders."
                            )
                    if len(filenames):
                        datasets = [
                            dataset for dataset in 
                            config["dependencies"].get(
                                "{}.{}".format(relpaths.figures / label, ext), []
                            ) if type(dataset) is str and 
                            dataset.endswith(".zenodo")
                        ]
                        if label in figures:
                            for fn in filenames:
                                if fn not in figures[label]["files"]:
                                    figures[label]["files"].append(fn)
                            for ds in datasets:
                                if ds not in figures[label]["datasets"]:
                                    figures[label]["datasets"].append(ds)
                        else:
                            figures[label] = {
                                "script": script, 
                                "files": filenames, 
                                "datasets": datasets
                            }

        # Parse free-floating graphics
        filenames = []
        for graphic in root.findall("GRAPHICS"):
            src = relpaths.src.absolute()
            path = (src / graphicspath / graphic.text).resolve().relative_to(src)
            if path.parts[0] == "figures":
                filenames.append(str(relpaths.src / path))
            elif path.parts[0] == "static":
                continue
            elif path.name in files.special_figures:
                continue
            else:
                warnings.warn(
                    f"Figure `{graphic.text}` must be in either the `src/figures` or `src/static` folders."
                )
        if len(filenames):
            figures["unknown"] = {
                "script": files.unknown,
                "files": filenames,
                "datasets": []
            }

        # Store as JSON
        scripts = {"figures": figures}
        with open(relpaths.temp / "scripts.json", "w") as f:
            print(json.dumps(scripts, indent=4), file=f)
