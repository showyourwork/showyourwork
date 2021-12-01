rule download_auto:
    """
    Download a showyourwork-managed figure dependency from Zenodo.

    """
    message:
        "Downloading dependency file {output[0]} from Zenodo..."
    output:
        temp("{dependency}") if config["CI"] else "{dependency}",
        "{dependency}.zenodo"
    wildcard_constraints:
        dependency="{}".format("|".join(files.zenodo_files_auto))
    params:
        zenodo_url=lambda w: zenodo.zenodo_url[w.dependency],
        deposit_id=lambda w: zenodo.deposit_id[w.dependency],
        file_name=lambda w: zenodo.file_name[w.dependency],
    shell:
        " && ".join(
            [
                "REDIRECT_URL=$(curl -Ls -o /dev/null -w %{{url_effective}} https://{params.zenodo_url}/record/{params.deposit_id})",
                "VERSION_ID=${{REDIRECT_URL#*record/}}",
                "curl https://{params.zenodo_url}/record/${{VERSION_ID}}/files/{params.file_name} --output {output[0]}", 
                "echo 'https://{params.zenodo_url}/record/{params.deposit_id}' > {output[1]}"
            ]
        )