
figures = config["tree"]["figures"]
for figure_name in figures:

    figscript = figures[figure_name]["script"]
    graphics = figures[figure_name]["graphics"]
    datasets = figures[figure_name]["datasets"]
    dependencies = figures[figure_name]["dependencies"]
    command = figures[figure_name]["command"]

    if command is None:
        continue

    if figscript is None:
        figscript = []

    rule:
        """
        Generate a figure given a figure script and optional dependencies.

        This is the workhorse of ``showyourwork``, which generates all the
        figures in the article.

        """
        input:
            figscript,
            datasets,
            dependencies,
            "environment.yml",
        output:
            report(graphics, category="Figure")
        conda:
            (paths.user / "environment.yml").as_posix()
        params:
            command=command
        shell:
            "{params.command}"