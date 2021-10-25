rule upload:
    """
    Upload a figure dependency to Zenodo.
    
    """
    message:
        "Uploading dependency file {input[0]} to Zenodo..."
    input:
        "{dependency}"
    output:
        "{dependency}.zenodo"
    wildcard_constraints:
        dependency="{}".format("|".join(files.zenodo_files_auto))
    conda:
        posix(abspaths.user / "environment.yml")
    params:
        action="upload",
        file_name=lambda w: zenodo.file_name[w.dependency],
        file_path=lambda w: zenodo.file_path[w.dependency],
        deposit_title=lambda w: zenodo.deposit_title[w.dependency],
        deposit_description=lambda w: zenodo.deposit_description[w.dependency],
        deposit_creators=lambda w: zenodo.deposit_creators[w.dependency],
        sandbox=lambda w: zenodo.sandbox[w.dependency],
        token_name=lambda w: zenodo.token_name[w.dependency],
        script=lambda w: zenodo.script[w.dependency],
        repo_url="{}/tree/{}".format(get_repo_url(), get_repo_sha())
    script:
        "../scripts/zenodo.py"