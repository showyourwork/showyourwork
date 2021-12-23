rule dotgraph:
    """
    Generate a .gv graph of the build process in the DOT language.

    """
    message:
        "Generating a DOT graph of the build process..."
    input:
        figures
    output:
        temp("src/figures/dag.gv")
    params:
        FIGURES=relpaths.figures
    conda:
        "../envs/dag.yml"
    script:
        "../scripts/dag.py"


for ext in config["figexts"]:
    if ext.lower() != "png":
        rule:
            """
            Convert a figure file to a PNG thumbnail.

            """
            message:
                "Generating PNG thumnail for figure `{input}`..."
            input:
                "{figure}." + ext
            output:
                temp("{figure}.png")
            conda:
                "../envs/dag.yml"
            shell:
                "convert -resize 500x500 {input} {output}"


rule dag:
    """
    Generate a DAG (directed acyclic graph) of the build process.

    """
    message:
        "Generating the DAG of the build process..."
    input:
        "src/figures/dag.gv",
        lambda w: [fig[: -len(Path(fig).suffix)] + ".png" for fig in figures(w) if Path(fig).name != "dag.pdf"]
    output:
        "src/figures/dag.pdf"
    conda:
        "../envs/dag.yml"
    shell:
        "dot -Tpdf {input[0]} > {output}"