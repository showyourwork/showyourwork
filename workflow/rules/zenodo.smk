from utils import zenodo


xnum = 1
dnum = 1
for host in ["zenodo", "zenodo_sandbox"]:


    # URL for the host
    zenodo_url = zenodo.zenodo_url[host]


    for deposit_id, entry in config[host].items():


        # Get entry metadata
        contents = entry["contents"]
        zip_files = entry["zip_files"]


        # Rules to download files individually
        for remote_file, local_file in contents.items():
            rulename = f"syw__download{dnum}"
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


        # Rules for zip files and tarballs
        for zip_file, zip_contents in zip_files.items():


            # Path to the local version of the zip file
            local_zip_file = contents[zip_file]


            # Rules to extract files individually
            for compressed_file, extracted_file in zip_contents.items():
                rulename = f"syw__extract{xnum}"
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