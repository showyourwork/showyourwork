import os
import json


localrules: metadata


rule metadata:
    """
    Generates article metadata from the output of the ``repo``
    and ``script_info`` rules. Saves it as a JSON in the temporary
    ``showyourwork`` directory.

    """
    message:
        "Generating article metadata..."
    input:
        [
            posix(file)
            for file in (relpaths.dot_github / "workflows").glob("showyourwork.yml")
        ],
        posix(relpaths.temp / "repo.json"),
        posix(relpaths.temp / "scripts.json"),
        dot_zenodo_files
    output:
        posix(relpaths.temp / "meta.json"),
    run:
        # Load the metadata
        with open(relpaths.temp / "repo.json", "r") as f:
            repo = json.load(f)
        with open(relpaths.temp / "scripts.json", "r") as f:
            scripts = json.load(f)
        meta = dict(repo=repo)

        # TODO: expose config by sharing it in meta, or construct this here?
        meta["repo"]["sha_tag_header"] = ""
        if config["style"]["show_git_sha_or_tag"]:
            if meta["repo"]["tag"] != "":
                meta["repo"]["sha_tag_header"] = (
                    f'Git tag: {meta["repo"]["tag"]}'
                )
            else:
                # The git default short hash is the first 7 characters:
                meta["repo"]["sha_tag_header"] = (
                    f'Git commit: {meta["repo"]["sha"][:7]}'
                )

        # Miscellaneous
        meta["CI"] = os.getenv("CI", "false") == "true"
        meta["run_id"] = os.getenv("GITHUB_RUN_ID", "")

        # Figure metadata
        meta["status"] = ScriptUpToDate
        meta["labels"] = {}
        for label, value in scripts["figures"].items():
            script = value["script"]
            if script != files.unknown:
                status = get_script_status(script)
                datasets = value["datasets"]
                meta["labels"]["{}_script".format(label)] = script
                meta["labels"]["{}_status".format(label)] = str(status)
                numbers = ["One", "Two", "Three"] # Built-in max of 3
                urls = []
                for dataset in datasets:
                    with open(dataset, "r") as f:
                        url = f.readlines()[0].replace("\n", "")
                        if url not in urls:
                            urls.append(url)
                for url, number in zip(urls, numbers):
                    meta["labels"]["{}_dataset{}".format(label, number)] = url
                meta["status"] = max(meta["status"], status)
        meta["status"] = str(meta["status"])

        # Store as JSON
        with open(relpaths.temp / "meta.json", "w") as f:
            print(json.dumps(meta, indent=4), file=f)
