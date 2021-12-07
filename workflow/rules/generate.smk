rule generate:
    """
    Generate a figure dependency.
    
    """
    message:
        "Generating figure dependency file {output[0]}..."
    input:
        script_dependencies
    output:
        "{dependency}"
    wildcard_constraints:
        dependency="{}".format("|".join(files.zenodo_files_auto)),
    conda:
        posix(abspaths.user / "environment.yml")
    params:
        shell_cmd=shell_cmd
    shell:
        "{params.shell_cmd}"