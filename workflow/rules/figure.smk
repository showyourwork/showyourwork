rule figure:
    """
    Generate a figure given a figure script and optional dependencies.

    """
    message:
        "Generating figure `{output}`..."
    input:
        figure_script,
        figure_script_dependencies,
        "environment.yml",
    output:
        report("{figure}", category="Figure")
    wildcard_constraints:
        figure="src/figures/(.*?)\.{}".format("|".join(config["figexts"])),
    params:
        script_name=script_name,
        FIGURES=relpaths.figures,
        TEMP=relpaths.temp
    conda:
        posix(abspaths.user / "environment.yml")
    script:
        "../scripts/figure.py"
