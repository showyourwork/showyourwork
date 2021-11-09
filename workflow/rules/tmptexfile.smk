localrules: tmptexfile


rule tmptexfile:
    """
    Copies the manuscript to a temporary TeX file for the XML build step.

    """
    message:
        "Writing temporary tex file..."
    input:
        posix(relpaths.tex / config["ms"]),
    output:
        temp(posix(relpaths.tex / "{}.tex".format(files.tmp_xml))),
    run:
        with open(relpaths.tex / config["ms"], "r") as f:
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
            raise ValueError(r"Missing `\documentclass` in file `{}`.".format(config["ms"]))
        with open(relpaths.tex / "{}.tex".format(files.tmp_xml), "w") as f:
            f.writelines(lines)