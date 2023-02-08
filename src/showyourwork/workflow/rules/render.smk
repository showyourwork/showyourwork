import json

from showyourwork import paths

SUPPORTED_FIGURE_EXTENSIONS = [".pdf", ".png", ".jpg", ".jpeg", ".eps", ".svg"]

rule syw__render_dag_config:
    input:
        rules.syw__dag.output,
        ensure_manuscript_dependencies,
    output:
        temp(paths.work(config).root / "_render_dag_config.json")
    run:
        with open(output[0], "w") as f:
            json.dump(config, f)

def generated_files(*_):
    deps = _get_manuscript_dependencies()
    return [d for d in deps if Path(d).suffix in SUPPORTED_FIGURE_EXTENSIONS]

rule syw__render_thumbnails:
    input:
        rules.syw__dag.output,
        files=generated_files,
        script=paths.package_data("showyourwork", "workflow", "scripts", "render_thumbnails.py"),
    output:
        directory(paths.work(config).root / "thumbnails")
    params:
        repo_path=paths.repo(config).root
    conda:
        paths.package_data("showyourwork", "workflow", "envs", "render_thumbnails.yml")
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
        script=paths.package_data("showyourwork", "workflow", "scripts", "render_dag.py")
    output:
        paths.work(config).subdir("dag") / "dag.pdf"
    params:
        repo_path=paths.repo(config).root,
        work_path=paths.work(config).root
    conda:
        paths.package_data("showyourwork", "workflow", "envs", "render_dag.yml")
    shell:
        """
        python {input.script:q} \\
            --config {input.config:q} \\
            --repo-path {params.repo_path:q} \\
            --work-path {params.work_path:q} \\
            --thumbnails-path {input.thumbnails:q} \\
            --output {output:q}
        """
