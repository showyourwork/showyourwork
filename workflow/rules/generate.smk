rule generate:
    """
    Generate a figure dependency.
    
    """
    message:
        "Generating figure dependency file {output[0]}..."
    input:
        lambda w: zenodo.generate_deps[w.dep_name]
    output:
        "src/figures/{dep_name}"
    wildcard_constraints:
        dep_name="{}".format("|".join(files.zenodo_files_auto)),
    conda:
        posix(abspaths.user / "environment.yml")
    params:
        generate_shell=lambda w: zenodo.generate_shell[w.dep_name]
    shell:
        "cd src/figures && {params.generate_shell}"