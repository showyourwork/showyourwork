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
        posix(relpaths.temp / "showyourwork.xml")
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
                graphicspath = Path(".")
            else:
                graphicspath = re.findall("\{(.*?)\}", graphicspath[-1].text)[0]
                graphicspath = Path(graphicspath)
        except:
            warnings.warn("Unable to parse the `graphicspath` command. Ignoring...")
            graphicspath = Path(".")

        # Parse labeled graphics inside `figure` environments
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
                        if str(path) == "figures/dag.pdf":
                            continue
                        elif path.parts[0] == "figures":
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

                        # Get all Zenodo datasets for this script
                        alldeps = config["dependencies"].get(
                            "{}.{}".format(relpaths.figures / label, ext), []
                        )
                        datasets = [
                            dataset for dataset in alldeps
                            if type(dataset) is str and
                            dataset.endswith(".zenodo")
                        ]

                        # Check if any of the datasets are inside a Zenodo
                        # tarball, and if so, add a link to that deposit
                        for tarball in zenodo.deposit_contents:
                            for f in zenodo.deposit_contents[tarball]:
                                if f in alldeps:
                                    if not f"{tarball}.zenodo" in datasets:
                                        datasets.append(f"{tarball}.zenodo")

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

        # Other figures
        other_filenames = set()

        # Parse graphics labeled with fig* inside `figure` environments
        # These won't be automatically generated, but we still
        # need to make them dependencies of the PDF
        for figure in root.findall("FIGURE"):
            labels = figure.findall("LABELSTAR")
            if len(labels) and labels[0].text is not None:
                for graphic in figure.findall("GRAPHICS"):
                    src = relpaths.src.absolute()
                    path = (src / graphicspath / graphic.text).resolve().relative_to(src)
                    if str(path) == "figures/dag.pdf":
                            continue
                    elif path.parts[0] == "figures":
                        other_filenames.add(str(relpaths.src / path))
                    elif path.parts[0] == "static":
                        continue
                    elif path.name in files.special_figures:
                        continue
                    else:
                        warnings.warn(
                            f"Figure `{graphic.text}` must be in either the `src/figures` or `src/static` folders."
                        )

        # Parse free-floating graphics
        for graphic in root.findall("GRAPHICS"):
            src = relpaths.src.absolute()
            path = (src / graphicspath / graphic.text).resolve().relative_to(src)
            if str(path) == "figures/dag.pdf":
                continue
            elif path.parts[0] == "figures":
                if path.parts[-1].startswith('orcid-id.'):
                    continue
                other_filenames.add(str(relpaths.src / path))
            elif path.parts[0] == "static":
                continue
            elif path.name in files.special_figures:
                continue
            else:
                warnings.warn(
                    f"Figure `{graphic.text}` must be in either the `src/figures` or `src/static` folders."
                )

        if len(other_filenames):
            figures["unknown"] = {
                "script": files.unknown,
                "files": list(other_filenames),
                "datasets": []
            }

        # Store as JSON
        scripts = {"figures": figures}
        with open(relpaths.temp / "scripts.json", "w") as f:
            print(json.dumps(scripts, indent=4), file=f)
