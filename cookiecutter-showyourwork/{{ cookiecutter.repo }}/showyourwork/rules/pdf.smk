localrules:
    texfile,
    stylesheet,


rule texfile:
    message:
        "Writing temporary tex file..."
    input:
        POSIX(TEX / "ms.tex"),
    output:
        temp(POSIX(TEX / "{}.tex".format(SYWTEXFILE))),
    run:
        with open(TEX / "ms.tex", "r") as f:
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
        with open(TEX / "{}.tex".format(SYWTEXFILE), "w") as f:
            f.writelines(lines)


rule stylesheet:
    message:
        "Generating stylesheet..."
    input:
        POSIX(TEMP / "meta.json"),
    output:
        temp(POSIX(TEX / "showyourwork.sty")),
    params:
        WORKFLOW=WORKFLOW,
        TEMP=TEMP,
        TEX=TEX,
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/stylesheet.py"


rule pdf:
    message:
        "Building pdf..."
    input:
        POSIX(TEX / "{}.tex".format(SYWTEXFILE)),
        [POSIX(file) for file in TEX.glob("*.bib")],
        POSIX(TEMP / "meta.json"),
        POSIX(TEX / "showyourwork.sty"),
        AUXFILES,
        class_files,
        figures,
    output:
        "ms.pdf",
    params:
        verbose=verbose,
        TEMP=TEMP,
        TEX=TEX,
        SYWTEXFILE=SYWTEXFILE,
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/pdf.py"
