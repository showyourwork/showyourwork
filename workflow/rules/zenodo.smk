znum = 1
for host in ["zenodo", "zenodo_sandbox"]:

    zenodo_url = zenodo.zenodo_url[host]

    for deposit_id, entry in config[host].items():

        # Get entry metadata
        token_name = entry["token_name"]
        id_type = entry["id_type"]
        contents = entry["contents"]
        zip_files = entry["zip_files"]

        # Rules for downloading the deposit contents
        if id_type == "version":

            # Loop over all loose files in the deposit
            for remote_file, local_file in contents.items():

                # Rule to download the file
                rulename = f"zenodo{znum}"
                znum += 1
                rule:
                    """
                    Download a Zenodo-hosted file.

                    """
                    name:
                        rulename
                    message:
                        "Downloading {output} from Zenodo..."
                    output:
                        report(local_file, category="Dataset")
                    params:
                        zenodo_url=zenodo_url,
                        deposit_id=deposit_id,
                        remote_file=remote_file
                    shell:
                        "curl https://{params.zenodo_url}/record/{params.deposit_id}/files/{params.remote_file} --output {output[0]}"

            # Loop over files contained in zip files & tarballs
            for zip_file, zip_contents in zip_files.items():

                # Path to the local version of the zip file
                local_zip_file = contents[zip_file]

                for compressed_file, extracted_file in zip_contents.items():

                    # Rule to extract the file
                    rulename = f"zenodo{znum}"
                    znum += 1
                    rule:
                        """
                        Extract a dataset from a zip file or tarball downloaded from Zenodo.

                        """
                        name:
                            rulename
                        message:
                            "Extracting {output}..."
                        input:
                            local_zip_file
                        output:
                            report(extracted_file, category="Dataset")
                        params:
                            compressed_file=compressed_file
                        script:
                            "../scripts/extract.py"
                            


        # TODO: Implement logic to re-download a dataset if the local
        # one is from a different version (track it in the config.json file)
        # Concept ids only.

        # TODO: Implement logic to handle upload/download of concept id datasets