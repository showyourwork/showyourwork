from utils import zenodo


xnum = 1
dnum = 1
cnum = 1
unum = 1
for host in ["zenodo", "zenodo_sandbox"]:


    # URL for the host
    zenodo_url = zenodo.zenodo_url[host]


    for deposit_id, entry in config[host].items():


        # Get entry metadata
        token_name = entry["token_name"]
        id_type = entry["id_type"]
        owner_ids = entry["owner_ids"]
        user_id = entry["user_id"]
        contents = entry["contents"]
        zip_files = entry["zip_files"]


        # Keep track of rules for this entry
        # so we can use `ruleorder` on them later
        entry["rules"] = []


        # Rule to upload all files in the deposit
        rulename = f"upload{unum}"
        entry["rules"].append(rulename)
        unum += 1
        rule:
            """
            Upload a deposit to Zenodo.

            """
            name:
                rulename
            message:
                "Uploading files to Zenodo..."
            input:
                list(contents.values())
            output:
                temp(touch(entry["upload_complete"]))
            conda:
                "../envs/main.yml"
            script:
                "../scripts/upload.py"


        # Rules to download files individually
        for remote_file, local_file in contents.items():
            rulename = f"download{dnum}"
            entry["rules"].append(rulename)
            dnum += 1
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
                script:
                    "../scripts/download.py"


        # Rules for zip files and tarballs (TODO: Hashing logic)
        for zip_file, zip_contents in zip_files.items():


            # Path to the local version of the zip file
            local_zip_file = contents[zip_file]


            # Rule to compress the files
            rulename = f"compress{cnum}"
            entry["rules"].append(rulename)
            cnum += 1
            rule:
                """
                Compress files into a zip file or tarball for a Zenodo upload.

                """
                name:
                    rulename
                message:
                    "Compressing {output}..."
                input:
                    list(zip_contents.values())
                output:
                    local_zip_file
                params:
                    compressed_files=list(zip_contents.keys())
                conda:
                    "../envs/main.yml"
                script:
                    "../scripts/compress.py"


            # Rules to extract files individually
            for compressed_file, extracted_file in zip_contents.items():
                rulename = f"extract{xnum}"
                entry["rules"].append(rulename)
                xnum += 1
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
                    conda:
                        "../envs/main.yml"
                    script:
                        "../scripts/extract.py"