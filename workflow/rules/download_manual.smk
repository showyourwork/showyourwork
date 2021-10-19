rule download_manual:
    """
    Download a figure dependency that was manually uploaded to Zenodo.

    """
    message:
        "Downloading dependency file {output[0]} from Zenodo..."
    output:
        temp("src/figures/{dep_name}") if config["CI"] else "src/figures/{dep_name}",
        "src/figures/{dep_name}.zenodo"
    wildcard_constraints:
        dep_name="{}".format("|".join(files.zenodo_files_manual))
    params:
        zenodo_url=lambda w: zenodo.zenodo_url[w.dep_name],
        zenodo_id=lambda w: zenodo.zenodo_id[w.dep_name],
    shell:
        " && ".join(
            [
                "curl https://{params.zenodo_url}/record/{params.zenodo_id}/files/{wildcards.dep_name} --output {output[0]}", 
                "echo 'https://{params.zenodo_url}/record/{params.zenodo_id}' > {output[1]}"
            ]
        )