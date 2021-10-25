rule download_manual:
    """
    Download a figure dependency that was manually uploaded to Zenodo.

    """
    message:
        "Downloading dependency file {output[0]} from Zenodo..."
    output:
        temp("src/{dependency}") if config["CI"] else "src/{dependency}",
        "src/{dependency}.zenodo"
    wildcard_constraints:
        dependency="{}".format("|".join(files.zenodo_files_manual))
    params:
        zenodo_url=lambda w: zenodo.zenodo_url[w.dependency],
        zenodo_id=lambda w: zenodo.zenodo_id[w.dependency],
        file_name=lambda w: zenodo.file_name[w.dependency]
    shell:
        " && ".join(
            [
                "curl https://{params.zenodo_url}/record/{params.zenodo_id}/files/{params.file_name} --output {output[0]}", 
                "echo 'https://{params.zenodo_url}/record/{params.zenodo_id}' > {output[1]}"
            ]
        )