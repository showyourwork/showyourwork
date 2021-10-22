localrules: stylesheet


rule stylesheet:
    """
    Generates the ``showyourwork.sty`` LaTeX style sheet from a Jinja
    template and copies it to the ``src`` directory.

    """
    message:
        "Generating stylesheet..."
    input:
        posix(relpaths.temp / "meta.json"),
    output:
        temp(posix(relpaths.tex / "showyourwork.sty")),
    params:
        WORKFLOW=abspaths.workflow,
        TEMP=relpaths.temp,
        TEX=relpaths.tex
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/stylesheet.py"