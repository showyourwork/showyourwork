
include: ".helpers.smk"
include: "Snakefile"

rule pdf:
    input:
        "tex/ms.tex",
        glob("tex/*.bib"),
        TEMP / "meta.json",
        figures
    output:
        "ms.pdf"
    run:
        run_pdf()


rule figure:
    input:
        figure_script,
        figure_data,
        "environment.yml"
    output:
        "{figure}"
    wildcard_constraints:
        figure="figures/(.*?)\.{}".format("|".join(FIGURE_EXTENSIONS))
    params:
        script_name=script_name,
        cached_figure=cached_figure,
        cache_cmd=cache_cmd
    conda:
        "environment.yml"
    shell:
        ("cd figures && "
        "[ -f {params.cached_figure} ] && "
        "{{ mv {params.cached_figure} . ; }} || "
        "{{ python {params.script_name} && {params.cache_cmd} ; }}")


rule repo_info:
    output:
        TEMP / "repo.json"
    priority: 
        99
    run:
        run_repo_info()


checkpoint script_info:
    input:
        "tex/ms.tex"
    output:
        TEMP / "scripts.json"
    run:
        run_script_info()


rule metadata:
    input:
        TEMP / "repo.json",
        TEMP / "scripts.json"
    output:
        TEMP / "meta.json"
    run:
        run_metadata()