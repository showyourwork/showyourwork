dnum = 1
unum = 1
for host in ["zenodo", "zenodo_sandbox"]:

    zenodo_url = zenodo.zenodo_url[host]
    entry = config[host]

    for dataset in entry:

        # Get entry metadata
        deposit_id = entry[dataset]["id"]
        token_name = entry[dataset]["token_name"]
        id_type = entry[dataset]["id_type"]
        data_script = entry[dataset]["script"]
        if data_script is None:
            data_script = []
        file_name = entry[dataset]["file_name"]
        file_path = entry[dataset]["file_path"]
        title = entry[dataset]["title"]
        description = entry[dataset]["description"]
        creators = entry[dataset]["creators"]
        contents = entry[dataset]["contents"]

        # Rule for downloading the deposit
        if id_type == "version":

            # User-friendly rule name
            rulename = f"download{dnum}"
            dnum += 1

            rule:
                """
                Download a Zenodo-hosted dataset.

                """
                name:
                    rulename
                message:
                    "Downloading `{output}` from Zenodo..."
                output:
                    report(dataset, category="Dataset")
                params:
                    zenodo_url=zenodo_url,
                    deposit_id=deposit_id,
                    file_name=file_name
                shell:
                    "curl https://{params.zenodo_url}/record/{params.deposit_id}/files/{params.file_name} --output {output[0]}"

        # TODO: Implement logic to re-download a dataset if the local
        # one is from a different version (track it in the config.json file)

        # TODO: Implement logic to handle upload/download of concept id datasets

        # TODO: Implement logic to handle tarballs