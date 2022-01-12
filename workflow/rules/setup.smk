rule setup:
    """
    Workflow setup.

    """
    input:
        (paths.temp / "tree.json").as_posix()
    conda:
        "../envs/setup.yml"
    script:
        "../scripts/setup.py"