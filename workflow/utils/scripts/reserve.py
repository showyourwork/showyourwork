import requests
import os


def reserve():
    """
    Pre-reserve a concept DOI on Zenodo.

    This is a user-facing command line utility called from the Makefile.
    Not to be run on GitHub Actions!

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
    if sandbox == "y":
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
        access_token = os.getenv(token_name, None)
        if access_token is None:
            print(f"Access token {token_name} not found.")
        else:
            break

    # Create a new deposit
    print("Pre-reserving Zenodo concept DOI...")
    r = requests.post(
        f"https://{zenodo_url}/api/deposit/depositions",
        params={
            "access_token": access_token,
        },
        json={},
    )
    if r.status_code > 204:
        data = r.json()
        for error in data.get("errors", []):
            data["message"] += " " + error["message"]
        raise Exception("Zenodo error {}: {}".format(data["status"], data["message"]))

    if r is not None:

        # Get the deposit id
        DEPOSIT_ID = r.json()["conceptrecid"]

        # We're done
        print(f"Deposit reserved on {zenodo_url} with concept DOI id {DEPOSIT_ID}.")


if __name__ == "__main__":
    reserve()