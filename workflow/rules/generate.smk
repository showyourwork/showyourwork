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
        path=lambda w: str(Path(zenodo.script[w.dependency]).parent),
        script=lambda w: str(Path(zenodo.script[w.dependency]).name)
    shell:
        "cd {params.path} && python {params.script}"