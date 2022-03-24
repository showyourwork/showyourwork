import sys
import json
import re
from pathlib import Path
from collections.abc import MutableMapping
from xml.etree.ElementTree import parse as ParseXMLTree


# Snakemake config (available automagically)
config = snakemake.config  # type:ignore
if config["showyourwork_path"]:
    sys.path.insert(1, config["showyourwork_path"])


# Import showyourwork
from showyourwork import paths, exceptions, zenodo
from showyourwork.tex import compile_tex


def flatten_zenodo_contents(d, parent_key="", default_path=None):
    """
    Flatten the `contents` dictionary of a Zenodo entry, filling
    in default mappings and removing zipfile extensions from the target path.

    Adapted from https://stackoverflow.com/a/6027615

    """
    if not default_path:
        default_path = paths.user().data.relative_to(paths.user().repo)
    items = []
    for k, v in d.items():
        new_key = (Path(parent_key) / k).as_posix() if parent_key else k
        if isinstance(v, MutableMapping):
            items.extend(
                flatten_zenodo_contents(
                    v, new_key, default_path=default_path
                ).items()
            )
        else:
            if v is None:
                # Use the default path
                # If inside a zipfile, remove the zipfile extension
                # from the target path
                zip_file = Path(new_key).parts[0]
                for ext in zenodo.zip_exts:
                    if zip_file.endswith(f".{ext}"):
                        mod_key = Path(new_key).parts[0][
                            : -len(f".{ext}")
                        ] / Path(*Path(new_key).parts[1:])
                        v = (Path(default_path) / mod_key).as_posix()
                        break
                else:
                    v = (Path(default_path) / new_key).as_posix()
            items.append((new_key, v))
    return dict(items)


def parse_overleaf():
    # Make sure `id` is defined
    config["overleaf"]["id"] = config["overleaf"].get("id", None)

    # Make sure `auto-commit` is defined
    config["overleaf"]["auto-commit"] = config["overleaf"].get(
        "auto-commit", False
    )

    # Make sure `push` and `pull` are defined and they are lists
    config["overleaf"]["push"] = config["overleaf"].get("push", [])
    if config["overleaf"]["push"] is None:
        config["overleaf"]["push"] = []
    elif type(config["overleaf"]["push"]) is not list:
        # TODO
        raise exceptions.ConfigError()
    config["overleaf"]["pull"] = config["overleaf"].get("pull", [])
    if config["overleaf"]["pull"] is None:
        config["overleaf"]["pull"] = []
    elif type(config["overleaf"]["pull"]) is not list:
        # TODO
        raise exceptions.ConfigError()

    # Ensure all files in `push` and `pull` are in the `src/tex` directory
    for file in config["overleaf"]["push"] + config["overleaf"]["pull"]:
        if not Path(file).resolve().is_relative_to(paths.user().tex):
            # TODO
            raise exceptions.ConfigError()

    # Ensure no overlap between `push` and `pull`.
    # User could in principle provide a directory in one
    # and a file within that directory in the other and that would
    # not trigger this error; we'll just have to let them live
    # dangerously!
    push_files = set(
        [
            str(Path(file).resolve().relative_to(paths.user().tex))
            for file in config["overleaf"]["push"]
        ]
    )
    pull_files = set(
        [
            str(Path(file).resolve().relative_to(paths.user().tex))
            for file in config["overleaf"]["pull"]
        ]
    )
    if len(push_files & pull_files):
        # TODO
        raise exceptions.ConfigError()


def parse_zenodo_datasets():
    """
    Parse the `zenodo` and `zenodo_sandbox` keys in the config file and
    populate entries with custom metadata.

    """
    for host in ["zenodo", "zenodo_sandbox"]:

        if host == "zenodo":
            tmp_path = paths.user().zenodo.relative_to(paths.user().repo)
        else:
            tmp_path = paths.user().zenodo_sandbox.relative_to(
                paths.user().repo
            )

        for deposit_id, entry in config[host].items():

            try:
                deposit_id = int(deposit_id)
            except:
                raise exceptions.ZenodoRecordNotFound(deposit_id)

            # Require that this is a static *version* ID
            entry["id_type"] = zenodo.get_id_type(
                deposit_id=deposit_id, zenodo_url=zenodo.zenodo_url[host]
            )
            if entry["id_type"] != "version":
                # TODO
                raise exceptions.InvalidZenodoIdType()

            # Deposit contents
            entry["destination"] = entry.get(
                "destination",
                str(paths.user().data.relative_to(paths.user().repo)),
            )
            contents = flatten_zenodo_contents(
                entry.get("contents", {}), default_path=entry["destination"]
            )

            # Handle files inside zipfiles, tarballs, etc.
            entry["zip_files"] = {}
            for source in list(contents.keys()):

                # Ensure the target is not a list
                target = contents[source]
                if type(target) is list:
                    # TODO
                    raise exceptions.ZenodoContentsError()

                # If it's a zipfile, add it to a separate entry in the config
                zip_file = Path(source).parts[0]
                if any(
                    [zip_file.endswith(f".{ext}") for ext in zenodo.zip_exts]
                ):
                    new_source = Path(*Path(source).parts[1:]).as_posix()
                    if zip_file in entry["zip_files"].keys():
                        entry["zip_files"][zip_file].update(
                            {new_source: target}
                        )
                    else:
                        entry["zip_files"][zip_file] = {new_source: target}

                    # Remove it from the `contents` entry
                    del contents[source]

                    # We'll host the zipfile in a temporary directory
                    contents[zip_file] = (
                        tmp_path / str(deposit_id) / zip_file
                    ).as_posix()

            entry["contents"] = contents


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
            # TODO
            raise exceptions.FigureFormatError(
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
                # TODO
                raise exceptions.FigureFormatError(
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
                # TODO
                raise exceptions.FigureFormatError(
                    "Command \marginicon must always come before the figure label."
                )

    # Check that there is exactly one label
    if len(labels) >= 2:
        # TODO
        raise exceptions.FigureFormatError(
            "A figure has multiple labels: `{}`".format(
                ", ".join(label.text for label in labels)
            )
        )
    elif len(labels) == 0:
        if len(figure.findall("LABELSTAR")) == 0:
            # TODO
            raise exceptions.FigureFormatError(
                "There is a figure without a label."
            )


def get_xml_tree():
    """"""
    # Parameters
    xmlfile = paths.user().preprocess / "showyourwork.xml"

    # Build the paper to get the XML file
    compile_tex(
        config,
        args=[
            "-r",
            "0",
        ],
        output_dir=paths.user().preprocess,
        stylesheet=paths.showyourwork().resources
        / "styles"
        / "preprocess.tex",
    )

    # Add <HTML></HTML> tags to the XML file
    if xmlfile.exists():
        with open(xmlfile, "r") as f:
            contents = f.read()
    else:
        raise exceptions.MissingXMLFile(
            r"Article parsing failed. Did you forget to `\usepackage{showyourwork}`?"
        )

    contents = "<HTML>\n" + contents + "</HTML>"
    with open(xmlfile, "w") as f:
        print(contents, file=f)

    # Load the XML tree
    return ParseXMLTree(paths.user().preprocess / "showyourwork.xml").getroot()


def get_json_tree():
    """"""
    # Get the XML article tree
    xml_tree = get_xml_tree()

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
        raise exceptions.GraphicsPathError()

    # Parse labeled graphics inside `figure` environments
    figures = {}
    for figure in xml_tree.findall("FIGURE"):

        # Ensure the figure environment conforms to the standard
        check_figure_format(figure)

        # Find all graphics included in this figure environment
        graphics = [
            str(
                (paths.user().tex / graphicspath / graphic.text)
                .resolve()
                .relative_to(paths.user().repo)
            )
            for graphic in figure.findall("GRAPHICS")
        ]

        # Are these static figures?
        static = all(
            [
                (paths.user().repo / graphic).parents[0]
                == paths.user().figures
                and (paths.user().static_figures / Path(graphic).name).exists()
                for graphic in graphics
            ]
        )

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
            for ext in config["script_extensions"]:
                script = paths.user().figure_scripts / f"{label}.{ext}"
                if script.exists():
                    script = str(script.relative_to(paths.user().repo))
                    command = config["scripts"][ext]
                    break
            else:
                script = None
                if static:
                    srcs = " ".join(
                        [
                            str(
                                (
                                    paths.user().static_figures
                                    / Path(graphic).name
                                ).relative_to(paths.user().repo)
                            )
                            for graphic in graphics
                        ]
                    )
                    dest = paths.user().figures.relative_to(paths.user().repo)
                    command = f"cp {srcs} {dest}"
                else:
                    command = None

        elif len(label_stars) and label_stars[0].text is not None:

            # User defined a \label{fig*:XXX}
            label = labels[0].text

            # There is no script associated with this figure
            script = None

            # But we may need to copy it over from the static dir
            if static:
                srcs = " ".join(
                    [
                        str(
                            (
                                paths.user().static_figures
                                / Path(graphic).name
                            ).relative_to(paths.user().repo)
                        )
                        for graphic in graphics
                    ]
                )
                dest = paths.user().figures.relative_to(paths.user().repo)
                command = f"cp {srcs} {dest}"
            else:
                command = None

        else:

            raise exceptions.FigureFormatError(
                "There is a figure without a label."
            )

        # Collect user-defined dependencies
        dependencies = config["dependencies"].get(script, [])

        # If any of the dependencies exists in a Zenodo deposit, infer
        # its URL here so we can add margin links to the PDF
        datasets = []
        for host in ["zenodo", "zenodo_sandbox"]:
            for deposit_id in config[host]:
                url = f"https://{zenodo.zenodo_url[host]}/record/{deposit_id}"
                for dep in dependencies:
                    if dep in config[host][deposit_id]["contents"].values():
                        datasets.append(url)
                    else:
                        for zip_file in config[host][deposit_id]["zip_files"]:
                            if (
                                dep
                                in config[host][deposit_id]["zip_files"][
                                    zip_file
                                ].values()
                            ):
                                datasets.append(url)
                                break
        datasets = list(set(datasets))

        # Format the command by replacing placeholders
        if command is not None:
            command = command.format(
                script=script,
                output=graphics,
                datasets=datasets,
                dependencies=dependencies,
            )

        # Add an entry to the tree
        figures[label] = {
            "script": script,
            "graphics": graphics,
            "datasets": datasets,
            "dependencies": dependencies,
            "command": command,
        }

    # Parse free-floating graphics
    graphics = [
        str(
            (paths.user().tex / graphicspath / graphic.text)
            .resolve()
            .relative_to(paths.user().repo)
        )
        for graphic in xml_tree.findall("GRAPHICS")
    ]

    # Ignore graphics that are dependencies of the texfile (such as orcid-ID.png)
    graphics = [
        graphic
        for graphic in graphics
        if graphic not in config["tex_files_out"]
    ]

    # Add an entry to the tree
    figures["free-floating"] = {
        "script": None,
        "graphics": graphics,
        "datasets": [],
        "dependencies": [],
        "command": None,
    }

    # The full tree (someday we'll have equations in here, too)
    tree = {"figures": figures}

    return tree


# Parse the `zenodo` key in the config
parse_zenodo_datasets()


# Parse overleaf config
parse_overleaf()


# Get the article tree
config["tree"] = get_json_tree()


# Make all of the graphics dependencies of the article
config["dependencies"][config["ms_tex"]] = config["dependencies"].get(
    config["ms_tex"], []
)
for figure_name in config["tree"]["figures"]:
    graphics = config["tree"]["figures"][figure_name]["graphics"]
    config["dependencies"][config["ms_tex"]].extend(
        [Path(graphic).as_posix() for graphic in graphics]
    )


# Gather the figure script & dataset info so we can access it on the TeX side
config["labels"] = {}
for label, value in config["tree"]["figures"].items():
    script = value["script"]
    if script is not None:
        config["labels"][f"{label}_script"] = script
    datasets = value["datasets"]
    # Note: built-in max of 3 datasets will be displayed
    for dataset, number in zip(datasets, ["One", "Two", "Three"]):
        config["labels"][f"{label}_dataset{number}"] = dataset


# Save the config file
with open(config["config_json"], "w") as f:
    print(json.dumps(config, indent=4), file=f)