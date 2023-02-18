import json

SUPPORTED_FIGURE_EXTENSIONS = [".pdf", ".png", ".jpg", ".jpeg", ".eps", ".svg"]
dag_directory = SYW__WORK_PATHS / "dag"

rule syw__render_dag_config:
    input:
        rules.syw__dag.output,
        ensure_all_document_dependencies,
    output:
        dag_directory / "_render_dag_config.json"
    run:
        with open(output[0], "w") as f:
            json.dump(config, f)

def generated_files(*_):
    deps = []
    for doc in SYW__DOCUMENTS:
        deps.extend(get_document_dependencies(doc)())
    return [d for d in set(deps) if Path(d).suffix in SUPPORTED_FIGURE_EXTENSIONS]

rule syw__render_thumbnails:
    input:
        rules.syw__dag.output,
        files=generated_files,
        script=package_data("showyourwork", "workflow", "scripts", "render_thumbnails.py"),
    output:
        directory(SYW__WORK_PATHS.root / "thumbnails")
    params:
        repo_path=SYW__REPO_PATHS.root
    conda:
        package_data("showyourwork", "workflow", "envs", "render_thumbnails.yml")
    shell:
        """
        python {input.script:q} \\
            --repo-path {params.repo_path:q} \\
            --output {output:q} \\
            {input.files:q}
        """

rule syw__render_dag:
    input:
        config=rules.syw__render_dag_config.output[0],
        thumbnails=rules.syw__render_thumbnails.output[0],
        script=package_data("showyourwork", "workflow", "scripts", "render_dag.py")
    output:
        config.get("dag", {}).get("render_to", "dag.pdf")
    params:
        repo_path=SYW__REPO_PATHS.root,
        work_path=SYW__WORK_PATHS.root
    conda:
        package_data("showyourwork", "workflow", "envs", "render_dag.yml")
    shell:
        """
        python {input.script:q} \\
            --config {input.config:q} \\
            --repo-path {params.repo_path:q} \\
            --work-path {params.work_path:q} \\
            --thumbnails-path {input.thumbnails:q} \\
            --output {output:q}
        """
