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
        # Note that the `deposit_id` variable is always a *concept* id
        # for showyourwork-managed figure dependencies.
        # To download the files, we need to resolve it to the latest
        # version id, which we can do by using `curl` to tell us what
        # URL we get redirected to when we access the webpage for the
        # concept id (REDIRECT_URL); the latest version id is then just
        # the last bit of the url (VERSION_ID).
        # We could do this a bit more elegantly using the REST API, though.
        " && ".join(
            [
                "REDIRECT_URL=$(curl -Ls -o /dev/null -w %{{url_effective}} https://{params.zenodo_url}/record/{params.deposit_id})",
                "VERSION_ID=${{REDIRECT_URL#*record/}}",
                "curl https://{params.zenodo_url}/record/${{VERSION_ID}}/files/{params.file_name} --output {output[0]}", 
                "echo 'https://{params.zenodo_url}/record/{params.deposit_id}' > {output[1]}"
            ]
        )