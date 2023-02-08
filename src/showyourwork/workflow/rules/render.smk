import json

from showyourwork import paths

rule syw__render_dag_config:
    input:
        rules.syw__dag.output,
        ensure_manuscript_dependencies,
    output:
        temp(paths.work(config).root / "_render_dag_config.json")
    run:
        with open(output[0], "w") as f:
            json.dump(config, f)

rule syw__render_dag:
    input:
        config=rules.syw__render_dag_config.output[0],
        script=paths.package_data("showyourwork", "workflow", "scripts", "render_dag.py")
    output:
        "dag.pdf"
    params:
        repo_path=str(paths.repo(config).root),
        work_path=str(paths.work(config).root)
    conda:
        paths.package_data("showyourwork", "workflow", "envs", "render_dag.yml")
    shell:
        """
        python "{input.script}" \\
            --config "{input.config}" \\
            --repo-path "{params.repo_path}" \\
            --work-path "{params.work_path}" \\
            --output "{output}"
        """
