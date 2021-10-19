rule upload:
    """
    Upload a figure dependency to Zenodo.
    
    """
    message:
        "Uploading dependency file {input[0]} to Zenodo..."
    input:
        "src/figures/{dep_name}"
    output:
        "src/figures/{dep_name}.zenodo"
    wildcard_constraints:
        dep_name="{}".format("|".join(files.zenodo_files_auto))
    conda:
        posix(abspaths.user / "environment.yml")
    params:
        action="upload",
        file_name=lambda w: zenodo.file_name[w.dep_name],
        file_path=lambda w: zenodo.file_path[w.dep_name],
        deposit_title=lambda w: zenodo.deposit_title[w.dep_name],
        deposit_description=lambda w: zenodo.deposit_description[w.dep_name],
        deposit_creators=lambda w: zenodo.deposit_creators[w.dep_name],
        sandbox=lambda w: zenodo.sandbox[w.dep_name],
        token_name=lambda w: zenodo.token_name[w.dep_name],
        generate_shell=lambda w: zenodo.generate_shell[w.dep_name],
        repo_url="{}/tree/{}".format(get_repo_url(), get_repo_sha())
    script:
        "../scripts/zenodo.py"