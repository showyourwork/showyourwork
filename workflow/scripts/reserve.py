"""
Pre-reserve a Zenodo DOI.

"""
import requests
import os


def check_status(r):
    if r.status_code > 204:
        data = r.json()
        for error in data.get("errors", []):
            data["message"] += " " + error["message"]
        print("Zenodo error {}: {}".format(data["status"], data["message"]))
        return None
    else:
        return r


def reserve(sandbox=False, token_name="ZENODO_TOKEN"):

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

    # Zenodo access token
    while True:
        try:
            token_name = input(
                "Name of environment variable containing Zenodo API token [ZENODO_TOKEN]: "
            )
        except (EOFError, KeyboardInterrupt):
            print("[Interrupted.]")
            return
        if token_name == "":
            token_name = "ZENODO_TOKEN"
        access_token = os.getenv(token_name, None)
        if access_token is None:
            print(
                f"Zenodo access token `{token_name}` not found. This should be set as an environment variable and/or repository secret."
            )
        else:
            break

    # Create a new deposit
    print("Pre-reserving Zenodo DOI...")
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
        DEPOSIT_ID = r.json()["id"]

        # We're done
        print(f"Deposit reserved on {zenodo_url} with id {DEPOSIT_ID}.")


if __name__ == "__main__":

    reserve()