
include: ".helpers.smk"
include: "Snakefile"


rule pdf:
    message:
        "Building pdf..."
    input:
        "tex/ms.tex",
        glob("tex/*.bib"),
        TEMP / "meta.json",
        ms_dependencies,
        figures
    output:
        "ms.pdf"
    run:
        run_pdf()


rule figure:
    message:
        "Generating figure `{output}`..."
    input:
        figure_script,
        figure_dependencies,
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


rule figure_extras:
    message:
        "Generating metadata for figure `{output}`..."
    output:
        temp(touch("{figure}.extras"))
    wildcard_constraints:
        figure="figures/(.*?)\.{}".format("|".join(FIGURE_EXTENSIONS))


rule repo_info:
    message:
        "Generating repo metadata..."
    output:
        TEMP / "repo.json"
    priority:
        99
    run:
        run_repo_info()


checkpoint script_info:
    message:
        "Building figure dependency graph..."
    input:
        "tex/ms.tex"
    output:
        TEMP / "scripts.json"
    run:
        run_script_info()


rule metadata:
    message:
        "Generating article metadata..."
    input:
        TEMP / "repo.json",
        TEMP / "scripts.json"
    output:
        TEMP / "meta.json"
    run:
        run_metadata()
