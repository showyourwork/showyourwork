rule download_ci:
    """
    Download a figure dependency from Zenodo that's managed by showyourwork.
    
    """
    message:
        "Downloading dependency file {output[0]} from Zenodo..."
    output:
        temp("src/figures/{dep_name}"),
        "src/figures/{dep_name}.zenodo"
    wildcard_constraints:
        dep_name="{}".format("|".join(files.zenodo_files_auto)),
    conda:
        posix(abspaths.user / "environment.yml")
    params:
        action="download",
        file_name=lambda w: zenodo.file_name[w.dep_name],
        file_path=lambda w: zenodo.file_path[w.dep_name],
        deposit_title=lambda w: zenodo.deposit_title[w.dep_name],
        sandbox=lambda w: zenodo.sandbox[w.dep_name],
        token_name=lambda w: zenodo.token_name[w.dep_name]
    script:
        "../scripts/zenodo.py"