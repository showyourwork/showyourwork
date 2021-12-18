rule dotgraph:
    """
    
    """
    input:
        figures
    output:
        temp("src/figures/dag.gv")
    conda:
        "../envs/dag.yml"
    script:
        "../scripts/dag.py"


for ext in config["figexts"]:
    if ext.lower() != "png":
        rule:
            """
    
            """
            input:
                "{figure}." + ext
            output:
                temp("{figure}.png")
            conda:
                "../envs/dag.yml"
            shell:
                "convert {input} {output}"


rule dag:
    """
    
    """
    input:
        "src/figures/dag.gv",
        lambda w: [fig[: -len(Path(fig).suffix)] + ".png" for fig in figures(w)]
    output:
        "src/figures/dag.pdf"
    conda:
        "../envs/dag.yml"
    shell:
        "dot -Tpdf {input[0]} > {output}"