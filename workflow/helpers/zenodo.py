"""
Main Zenodo interface.

"""
from .exceptions import ShowyourworkException
import requests
import os
import json
import socket
import requests.packages.urllib3.util.connection as urllib3_cn


# Force IPv4 for Zenodo connections. This is a temporary hack
# that sometimes greatly speeds up access to the API;
# TODO: we need investigate what's actually happening here!
# See https://stackoverflow.com/a/62599037
# and https://stackoverflow.com/a/46972341
urllib3_cn.allowed_gai_family = lambda: socket.AF_INET


def check_status(r):
    """
    Parse a requests return object and raise a custom exception
    for a >200-level status code.

    """
    if r.status_code > 204:
        data = r.json()
        for error in data.get("errors", []):
            data["message"] += " " + error["message"]
        raise ShowyourworkException(
            "Zenodo error {}: {}".format(data["status"], data["message"])
        )
    return r


def get_access_token(token_name):
    """
    Return the access token stored in the environment variable `token_name`.

    """
    access_token = os.getenv(token_name, None)
    if access_token is None:
        raise ValueError(
            f"Zenodo access token `{token_name}` not found. "
            "This should be set as an environment variable "
            "and/or a repository secret."
        )
    return access_token


def reserve():
    """
    Pre-reserve a concept DOI on Zenodo.
    This is a user-facing command line utility called from the Makefile.

    """
    # Zenodo Sandbox (for testing) or Zenodo?
    while True:
        try:
            sandbox = input("Pre-reserve on Zenodo Sandbox (y/n) [n]: ")
        except (EOFError, KeyboardInterrupt):
            print("[Interrupted.]")
            return
        if sandbox == "" or sandbox == "n":
            zenodo_url = "zenodo.org"
            break
        elif sandbox == "y":
            zenodo_url = "sandbox.zenodo.org"
            break
        else:
            print("Please enter either 'y' or 'n'.")

    # Default token name
    if sandbox:
        default_token_name = "ZENODO_SANDBOX_TOKEN"
    else:
        default_token_name = "ZENODO_TOKEN"

    # Zenodo access token
    while True:
        try:
            token_name = input(
                f"Name of environment variable containing Zenodo API token [{default_token_name}]: "
            )
        except (EOFError, KeyboardInterrupt):
            print("[Interrupted.]")
            return
        if token_name == "":
            token_name = default_token_name
        try:
            access_token = get_access_token(token_name)
        except Exception as e:
            print(str(e))
        else:
            break

    # Create a new deposit
    print("Pre-reserving Zenodo concept DOI...")
    r = check_status(
        requests.post(
            f"https://{zenodo_url}/api/deposit/depositions",
            params={
                "access_token": access_token,
            },
            json={},
        )
    )

    if r is not None:

        # Get the deposit id
        DEPOSIT_ID = r.json()["conceptrecid"]

        # We're done
        print(f"Deposit reserved on {zenodo_url} with concept DOI id {DEPOSIT_ID}.")


def resolve_id(
    deposit_id,
    zenodo_url="zenodo.org",
    token_name="ZENODO_TOKEN",
    max_pages=100,
    results_per_page=10,
):
    """
    Resolves a Zenodo (or Zenodo Sandbox) id into a concept id, a version
    id, and a draft id (if the user has the correct privileges).

    """

    # The dict we'll return below
    id_dict = {"input": int(deposit_id)}

    # Get the access token
    access_token = get_access_token(token_name)

    # Try to find a published record (no authentication needed)
    r = requests.get(f"https://{zenodo_url}/api/records/{deposit_id}")
    data = r.json()

    if r.status_code > 204:

        if "PID is not registered" in data.get("message", ""):

            try:

                # There's no public record; let's search for a draft
                # deposit with this concept or version ID (authentication needed)
                for page in range(1, max_pages + 1):
                    r = requests.get(
                        f"https://{zenodo_url}/api/deposit/depositions",
                        params={
                            "q": "",
                            "access_token": access_token,
                            "status": "draft",
                            "size": results_per_page,
                            "page": page,
                            "sort": "mostrecent",
                            "all_versions": 1,
                        },
                    )
                    data = r.json()
                    for deposit in data:
                        if int(deposit["conceptrecid"]) == deposit_id:

                            # We found a deposit draft with this concept id
                            id_dict["concept"] = int(deposit_id)
                            id_dict["draft"] = int(deposit["id"])
                            id_dict["version"] = None  # Doesn't exist yet
                            raise StopIteration

                        elif int(deposit["id"]) == deposit_id:

                            # We found a deposit draft with this version id
                            id_dict["concept"] = int(deposit["conceptrecid"])
                            id_dict["draft"] = int(deposit_id)
                            id_dict["version"] = None  # Doesn't exist yet
                            raise StopIteration

                    if len(data) == 0:

                        # We reached the end and found nothing
                        raise ShowyourworkException(
                            f"Record/deposit with id {deposit_id} not found.",
                            context=f"The provided `id` {deposit_id} does "
                            "not seem to be a valid Zenodo id.",
                            brief=f"Invalid Zenodo id {deposit_id}.",
                        )

            except StopIteration:

                # We found a record; we're good
                pass

        else:

            # Something unexpected happened...
            status = data.get("status", "unknown")
            message = data.get("message", "An error occurred while accessing Zenodo.")
            context = "Zenodo error {}: {}".format(status, message)
            raise ShowyourworkException(
                message,
                context=context,
                brief="An error occurred while accessing Zenodo.",
            )

    else:

        # This is either a version id or a concept id
        id_dict["concept"] = int(data["conceptrecid"])
        id_dict["version"] = int(data["id"])

        # Now if we have the correct authentication, let's get the
        # id of the latest draft (and create one if necessary).
        try:

            r = check_status(
                requests.get(
                    f"https://{zenodo_url}/api/deposit/depositions/{data['id']}",
                    params={"access_token": access_token},
                )
            )
            data = r.json()
            latest_draft = data["links"].get("latest_draft", None)
            if latest_draft:

                # A draft exists; get its id
                id_dict["draft"] = int(latest_draft.split("/")[-1])

            else:

                # There's no draft, so we need to create one
                r = check_status(
                    requests.post(
                        f"https://{zenodo_url}/api/deposit/depositions/{deposit_id}/actions/newversion",
                        params={"access_token": access_token},
                    )
                )
                id_dict["draft"] = int(r.json()["links"]["latest_draft"].split("/")[-1])

        except ShowyourworkException as e:

            if "401" in str(e) or "403" in str(e):
                # We don't have the permissions for this; no biggie
                id_dict["draft"] = None
            else:
                raise e

    return id_dict


def get_id_type(
    deposit_id,
    zenodo_url="zenodo.org",
    token_name="ZENODO_TOKEN",
    max_pages=100,
    results_per_page=10,
):
    """
    Determines whether a given Zenodo `id` corresponds to
    a concept id or a version id (and raises an error otherwise).

    """
    # Get the access token
    access_token = get_access_token(token_name)

    # Try to find a published record (no authentication needed)
    r = requests.get(f"https://{zenodo_url}/api/records/{deposit_id}")
    data = r.json()

    if r.status_code > 204:

        if "PID is not registered" in data.get("message", ""):

            # There's no public record; let's search for a draft
            # deposit with this concept or version ID (authentication
            # needed)
            for page in range(1, max_pages + 1):
                r = requests.get(
                    f"https://{zenodo_url}/api/deposit/depositions",
                    params={
                        "q": "",
                        "access_token": access_token,
                        "status": "draft",
                        "size": results_per_page,
                        "page": page,
                        "sort": "mostrecent",
                        "all_versions": 1,
                    },
                )
                data = r.json()
                for deposit in data:
                    if int(deposit["conceptrecid"]) == deposit_id:

                        # We found a deposit draft with this concept id
                        return "concept"

                    elif int(deposit["id"]) == deposit_id:

                        # We found a deposit draft with this version id
                        return "version"

                if len(data) == 0:

                    # We reached the end and found nothing
                    raise ShowyourworkException(
                        f"Record/deposit with id {deposit_id} not found.",
                        context=f"The provided `id` {deposit_id} does "
                        "not seem to be a valid concept or version Zenodo id.",
                        brief=f"Invalid Zenodo id {deposit_id}.",
                    )

        else:

            # Something unexpected happened...
            status = data.get("status", "unknown")
            message = data.get("message", "An error occurred while accessing Zenodo.")
            context = "Zenodo error {}: {}".format(status, message)
            raise ShowyourworkException(
                message,
                context=context,
                brief="An error occurred while accessing Zenodo.",
            )

    else:

        # This is a public record
        if int(deposit_id) == int(data["conceptrecid"]):
            return "concept"
        elif int(deposit_id) == int(data["id"]):
            return "version"
        else:
            raise ShowyourworkException(
                f"Record/deposit with id {deposit_id} not found.",
                context=f"The provided `id` {deposit_id} does "
                "not seem to be a valid concept or version Zenodo id.",
                brief=f"Invalid Zenodo id {deposit_id}.",
            )


def upload_simulation(
    file_name,
    draft_id,
    deposit_title,
    deposit_description,
    deposit_creators,
    zenodo_url="zenodo.org",
    token_name="ZENODO_TOKEN",
    file_path=".",
    script="",
    repo_url="",
):
    # Retrieve the access token
    access_token = get_access_token(token_name)

    # Get the bucket url for this draft so we can upload files
    # Also grab the latest version doi if we need to fall back to it below
    r = check_status(
        requests.get(
            f"https://{zenodo_url}/api/deposit/depositions/{draft_id}",
            params={"access_token": access_token},
        )
    )
    bucket_url = r.json()["links"]["bucket"]
    latest_id = r.json()["links"]["latest_html"].split("/")[-1]

    # Get the ID of the previously uploaded file (if it exists),
    # then delete it so we can upload a new version.
    print("Deleting old file(s)...")
    r = check_status(
        requests.get(
            f"https://{zenodo_url}/api/deposit/depositions/{draft_id}/files",
            params={"access_token": access_token},
        )
    )
    for file in r.json():
        if file["filename"] == file_name:
            FILE_ID = file["id"]
            r = check_status(
                requests.delete(
                    f"https://{zenodo_url}/api/deposit/depositions/{draft_id}/files/{FILE_ID}",
                    params={"access_token": access_token},
                )
            )
            break

    # Upload the new version of the file
    print("Uploading the file...")
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
    if script == "unknown-script":
        description = f"{deposit_description}<br/><br/>Created using <a href='https://github.com/rodluger/showyourwork'>showyourwork</a> from <a href='{repo_url}'>this GitHub repo</a>.<br/>"
    else:
        description = f"{deposit_description}<br/><br/>Created using <a href='https://github.com/rodluger/showyourwork'>showyourwork</a> from <a href='{repo_url}'>this GitHub repo</a> using the following command:<br/><pre><code class='language-bash'>cd src/figures && python {script}</code></pre><br/>"
    data = {
        "metadata": {
            "title": deposit_title,
            "upload_type": "dataset",
            "description": description,
            "creators": [{"name": name} for name in deposit_creators],
        }
    }
    r = check_status(
        requests.put(
            f"https://{zenodo_url}/api/deposit/depositions/{draft_id}",
            params={"access_token": access_token},
            data=json.dumps(data),
            headers={"Content-Type": "application/json"},
        )
    )

    # DEBUG
    # TODO: Remove me
    breakpoint()

    # Publish the deposit
    print("Publishing the deposit...")
    try:
        r = check_status(
            requests.post(
                f"https://{zenodo_url}/api/deposit/depositions/{draft_id}/actions/publish",
                params={"access_token": access_token},
            )
        )
    except ShowyourworkException as e:
        if "New version's files must differ from all previous versions" in str(e):
            print("No change in the deposit's files. Aborting.")
            # Revert to the previous deposit (the latest version) ID
            draft_id = latest_id
        else:
            raise e

    # Store the deposit URL
    deposit_url = f"https://{zenodo_url}/record/{draft_id}"
    with open(f"{file_path}/{file_name}.zenodo", "w") as f:
        print(deposit_url, file=f)