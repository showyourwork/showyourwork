"""
Upload or downlod a dataset to/from Zenodo.
This script is called from the ``zenodo`` rule.

"""
import requests
import os
import json


class ZenodoError(Exception):
    pass


def check_status(r):
    if r.status_code > 204:
        data = r.json()
        for error in data.get("errors", []):
            data["message"] += " " + error["message"]
        raise ZenodoError("Zenodo error {}: {}".format(data["status"], data["message"]))
    return r


def find_deposit(deposit_title, sandbox=False, token_name="ZENODO_TOKEN"):

    # Zenodo Sandbox (for testing) or Zenodo?
    if sandbox:
        zenodo_url = "sandbox.zenodo.org"
    else:
        zenodo_url = "zenodo.org"

    # Retrieve the access token
    access_token = os.getenv(token_name, None)
    if access_token is None:
        raise ValueError(
            f"Zenodo access token `{token_name}` not found. This should be set as an environment variable and/or repository secret."
        )

    # Search for an existing deposit with the given title
    print("Searching for existing deposit...")
    r = check_status(
        requests.get(
            f"https://{zenodo_url}/api/deposit/depositions",
            params={
                "q": deposit_title,
                "access_token": access_token,
            },
        )
    )
    deposit = None
    for entry in r.json():
        if entry["title"] == deposit_title:
            deposit = entry
            break

    return deposit


def upload_simulation(
    file_name,
    deposit_title,
    deposit_description,
    deposit_creators,
    sandbox=False,
    token_name="ZENODO_TOKEN",
    file_path=".",
    generate_shell="",
    repo_url="",
):

    # Upload to sandbox (for testing) or to actual Zenodo?
    if sandbox:
        zenodo_url = "sandbox.zenodo.org"
    else:
        zenodo_url = "zenodo.org"

    # Retrieve the access token
    access_token = os.getenv(token_name, None)
    if access_token is None:
        raise ValueError(
            f"Zenodo access token `{token_name}` not found. This should be set as an environment variable and/or repository secret."
        )

    # Search for an existing deposit with the given title
    deposit = find_deposit(deposit_title, sandbox=sandbox, token_name=token_name)

    # Either retrieve the deposit or create a new one
    if deposit:

        # Get the deposit id
        DEPOSIT_ID = deposit["id"]

        # Update the existing deposit
        print("Retrieving existing deposit...")
        try:

            # Create a new version draft
            r = check_status(
                requests.post(
                    f"https://{zenodo_url}/api/deposit/depositions/{DEPOSIT_ID}/actions/newversion",
                    params={"access_token": access_token},
                )
            )
            DEPOSIT_ID = r.json()["links"]["latest_draft"].split("/")[-1]

        except ZenodoError as e:

            if "403: Invalid action" in str(e):

                # Seems like we already have a draft. Let's use it
                DEPOSIT_ID = deposit["links"]["latest_draft"].split("/")[-1]

            else:

                raise e

        # Get the ID of the previously uploaded file (if it exists),
        # then delete it so we can upload a new version.
        print("Deleting old file...")
        r = check_status(
            requests.get(
                f"https://{zenodo_url}/api/deposit/depositions/{DEPOSIT_ID}/files",
                params={"access_token": access_token},
            )
        )
        for file in r.json():
            if file["filename"] == file_name:
                FILE_ID = file["id"]
                r = check_status(
                    requests.delete(
                        f"https://{zenodo_url}/api/deposit/depositions/{DEPOSIT_ID}/files/{FILE_ID}",
                        params={"access_token": access_token},
                    )
                )
                break

        # Upload the new version of the file
        print("Uploading new file...")
        with open(os.path.join(file_path, file_name), "rb") as fp:
            r = check_status(
                requests.post(
                    f"https://{zenodo_url}/api/deposit/depositions/{DEPOSIT_ID}/files",
                    data={"name": file_name},
                    files={"file": fp},
                    params={"access_token": access_token},
                )
            )

        # Add some metadata
        print("Adding metadata...")
        data = {
            "metadata": {
                "title": deposit_title,
                "upload_type": "dataset",
                "description": f"{deposit_description}<br/><br/>Created using <a href='https://github.com/rodluger/showyourwork'>showyourwork</a> from <a href='{repo_url}'>this GitHub repo</a> using the following command: <pre><code class='language-bash'>cd src/figures && {generate_shell}</code></pre>",
                "creators": [{"name": name} for name in deposit_creators],
            }
        }
        r = check_status(
            requests.put(
                f"https://{zenodo_url}/api/deposit/depositions/{DEPOSIT_ID}",
                params={"access_token": access_token},
                data=json.dumps(data),
                headers={"Content-Type": "application/json"},
            )
        )

        # Publish the deposit
        print("Publishing the deposit...")
        try:
            r = check_status(
                requests.post(
                    f"https://{zenodo_url}/api/deposit/depositions/{DEPOSIT_ID}/actions/publish",
                    params={"access_token": access_token},
                )
            )
        except ZenodoError as e:
            if "New version's files must differ from all previous versions" in str(e):
                print("No change in the deposit's files. Aborting.")
                # Revert to the previous deposit ID
                DEPOSIT_ID = deposit["links"]["latest_html"].split("/")[-1]
                pass
            else:
                raise e

    else:

        # Create a new deposit
        print("Creating a new deposit...")
        r = check_status(
            requests.post(
                f"https://{zenodo_url}/api/deposit/depositions",
                params={"access_token": access_token},
                json={},
            )
        )

        # Get the deposit id
        DEPOSIT_ID = r.json()["id"]

        # Upload the file
        print("Uploading the file...")
        bucket_url = r.json()["links"]["bucket"]
        with open(os.path.join(file_path, file_name), "rb") as fp:
            r = check_status(
                requests.put(
                    "%s/%s" % (bucket_url, file_name),
                    data=fp,
                    params={"access_token": access_token},
                )
            )

        # Add some metadata
        print("Adding metadata...")
        data = {
            "metadata": {
                "title": deposit_title,
                "upload_type": "dataset",
                "description": f"{deposit_description}<br/><br/>Created using <a href='https://github.com/rodluger/showyourwork'>showyourwork</a> from <a href='{repo_url}'>this GitHub repo</a> using the following command: <pre><code class='language-bash'>cd src/figures && {generate_shell}</code></pre>",
                "creators": [{"name": name} for name in deposit_creators],
            }
        }
        r = check_status(
            requests.put(
                f"https://{zenodo_url}/api/deposit/depositions/{DEPOSIT_ID}",
                params={"access_token": access_token},
                data=json.dumps(data),
                headers={"Content-Type": "application/json"},
            )
        )

        # Publish the deposit
        print("Publishing the deposit...")
        r = check_status(
            requests.post(
                f"https://{zenodo_url}/api/deposit/depositions/{DEPOSIT_ID}/actions/publish",
                params={"access_token": access_token},
            )
        )

    # Store the deposit URL
    deposit_url = f"https://{zenodo_url}/record/{DEPOSIT_ID}"
    with open(f"{file_path}/{file_name}.zenodo", "w") as f:
        print(deposit_url, file=f)


def download_simulation(
    file_name,
    deposit_title,
    sandbox=False,
    token_name="ZENODO_TOKEN",
    file_path=".",
):

    # Donwload from sandbox (for testing) or from actual Zenodo?
    if sandbox:
        zenodo_url = "sandbox.zenodo.org"
    else:
        zenodo_url = "zenodo.org"

    # Retrieve the access token
    access_token = os.getenv(token_name, None)
    if access_token is None:
        raise ValueError(
            f"Zenodo access token `{token_name}` not found. This should be set as an environment variable and/or repository secret."
        )

    # Search for an existing deposit with the given title
    deposit = find_deposit(deposit_title, sandbox=sandbox, token_name=token_name)
    if deposit is None:
        raise ZenodoError("Cannot find deposit with the given title.")

    # Get the latest *submitted* version
    if deposit["submitted"]:
        DEPOSIT_ID = deposit["id"]
    else:
        DEPOSIT_ID = deposit["links"]["latest_html"].split("/")[-1]

    # Download the file
    print("Downloading file...")
    r = check_status(
        requests.get(
            f"https://{zenodo_url}/api/deposit/depositions/{DEPOSIT_ID}/files",
            params={"access_token": access_token},
        )
    )
    for file in r.json():
        if file["filename"] == file_name:
            FILE_ID = file["id"]
            r = check_status(
                requests.get(
                    f"https://{zenodo_url}/api/deposit/depositions/{DEPOSIT_ID}/files/{FILE_ID}",
                    params={"access_token": access_token},
                )
            )
            url = r.json()["links"]["download"]
            r = check_status(
                requests.get(
                    url,
                    params={"access_token": access_token},
                )
            )
            with open(os.path.join(file_path, file_name), "wb") as f:
                f.write(r.content)
            break
    else:
        raise ZenodoError("Unable to download the file.")

    # Store the deposit URL
    deposit_url = f"https://{zenodo_url}/record/{DEPOSIT_ID}"
    with open(f"{file_path}/{file_name}.zenodo", "w") as f:
        print(deposit_url, file=f)


# Upload or download the file
if snakemake.params["action"] == "upload":
    upload_simulation(
        snakemake.params["file_name"],
        snakemake.params["deposit_title"],
        snakemake.params["deposit_description"],
        snakemake.params["deposit_creators"],
        sandbox=snakemake.params["sandbox"],
        token_name=snakemake.params["token_name"],
        file_path=snakemake.params["file_path"],
        repo_url=snakemake.params["repo_url"],
        generate_shell=snakemake.params["generate_shell"],
    )
else:
    download_simulation(
        snakemake.params["file_name"],
        snakemake.params["deposit_title"],
        sandbox=snakemake.params["sandbox"],
        token_name=snakemake.params["token_name"],
        file_path=snakemake.params["file_path"],
    )
