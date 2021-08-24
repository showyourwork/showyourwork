def check_figure_format(figure):
    """
    Check that all figures are declared correctly in `tex/ms.tex`
    so we can parse them corresponding XML tree.

    """
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


rule tmptexfile:
    message:
        "Writing temporary tex file..."
    input:
        posix(TEX / "ms.tex")
    output:
        temp(posix(TEX / "{}.tex".format(TMPTEXFILE)))
    run:
        with open(TEX / "ms.tex", "r") as f:
            lines = f.readlines()
        for idx, line in enumerate(lines):
            if line.startswith(r"\documentclass"):
                lines = (
                    lines[: idx + 1]
                    + [r"\usepackage{sywxml}" + "\n"]
                    + lines[idx + 1 :]
                )
                break
        else:
            raise ValueError(r"Missing `\documentclass` in file `tex/ms.tex`.")
        with open(TEX / "{}.tex".format(TMPTEXFILE), "w") as f:
            f.writelines(lines)


rule xml:
    message:
        "Generating XML tree..."
    input:
        class_files,
        aux_files,
        posix(TEX / "sywxml.sty"),
        posix(TEX / "{}.tex".format(TMPTEXFILE))
    output:
        temp(posix(TEMP / "showyourwork.xml")),
        temp(posix(TEMP / "{}.pdf".format(TMPTEXFILE)))
    run:
        tectonic_args = ["-r", "0", "-o", TEMP]
        if verbose:
            tectonic_args += ["--print"]
        else:
            tectonic_args += ["--chatter", "minimal"]
        subprocess.check_call(
            ["tectonic"] + tectonic_args + [TEX / "{}.tex".format(TMPTEXFILE)]
        )


checkpoint script_info:
    message:
        "Building figure dependency tree..."
    input:
        posix(TEMP / "showyourwork.xml")
    output:
        temp(posix(TEMP / "scripts.json"))
    run:
        # Load the XML tree
        root = ParseXMLTree(TEMP / "showyourwork.xml").getroot()

        # Parse graphics inside `figure` environments
        figures = {}
        for figure in root.findall("FIGURE"):
            check_figure_format(figure)
            labels = figure.findall("LABEL")
            if len(labels):
                label = labels[0].text
                if label is not None:
                    script = posix(FIGURES / "{}.py".format(label))
                    files = []
                    for graphic in figure.findall("GRAPHICS"):
                        if graphic.text.startswith("figures/"):
                            files.append(graphic.text.replace("/", os.sep))
                        else:
                            files.append(posix(FIGURES / graphic.text))
                    figures[label] = {"script": script, "files": files}

        # Parse free-floating graphics
        files = []
        for graphic in root.findall("GRAPHICS"):
            if graphic.text.startswith("figures/"):
                files.append(graphic.text.replace("/", os.sep))
            else:
                files.append(posix(FIGURES / graphic.text))
        figures["unknown"] = {
            "script": UNKNOWN_SCRIPT,
            "files": files,
        }

        # Store as JSON
        scripts = {"figures": figures}
        with open(TEMP / "scripts.json", "w") as f:
            print(json.dumps(scripts, indent=4), file=f)
