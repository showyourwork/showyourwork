rule download_ci:
    """
    Download a figure dependency from Zenodo that's managed by showyourwork.
    
    """
    message:
        "Downloading dependency file {output[0]} from Zenodo..."
    output:
        temp("{dependency}"),
        "{dependency}.zenodo"
    wildcard_constraints:
        dependency="{}".format("|".join(files.zenodo_files_auto)),
    conda:
        posix(abspaths.user / "environment.yml")
    params:
        action="download",
        file_name=lambda w: zenodo.file_name[w.dependency],
        file_path=lambda w: zenodo.file_path[w.dependency],
        deposit_title=lambda w: zenodo.deposit_title[w.dependency],
        sandbox=lambda w: zenodo.sandbox[w.dependency],
        token_name=lambda w: zenodo.token_name[w.dependency]
    script:
        "../scripts/zenodo.py"