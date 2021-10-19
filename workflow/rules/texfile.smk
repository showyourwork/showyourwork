localrules: texfile


rule texfile:
    """
    Copy the manuscript to a temporary TeX file for the final build.

    """
    message:
        "Writing temporary tex file..."
    input:
        posix(relpaths.tex / "ms.tex"),
    output:
        temp(posix(relpaths.tex / "{}.tex".format(files.tmp_syw))),
    run:
        with open(relpaths.tex / "ms.tex", "r") as f:
            lines = f.readlines()
        for idx, line in enumerate(lines):
            if line.startswith(r"\documentclass"):
                lines = (
                    lines[: idx + 1]
                    + [r"\usepackage{showyourwork}" + "\n"]
                    + lines[idx + 1 :]
                )
                break
        else:
            raise ValueError(r"Missing `\documentclass` in file `tex/ms.tex`.")
        with open(relpaths.tex / "{}.tex".format(files.tmp_syw), "w") as f:
            f.writelines(lines)