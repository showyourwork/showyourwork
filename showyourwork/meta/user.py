from ..constants import *
from ..utils import save_json
import questionary
import json


def get_user_metadata(clobber=True):

    if clobber or not (TEMP / "user" / "user.json").exists():

        # Custom style
        style = questionary.Style([])

        # Defaults
        if (TEMP / "user" / "user.json").exists():
            with open(TEMP / "user" / "user.json", "r") as f:
                user = json.load(f)
        else:
            user = {
                "name": "",
                "email": "",
                "affil": "",
                "orcid": "",
                "github": "",
            }

        # Command-line questionnaire
        user["name"] = questionary.text(
            "Full Name:",
            qmark="[ShowYourWork!]",
            default=user["name"],
            style=style,
        ).ask()

        user["email"] = questionary.text(
            "Email:",
            qmark="[ShowYourWork!]",
            default=user["email"],
            style=style,
        ).ask()

        user["affil"] = questionary.text(
            "Affiliation:",
            qmark="[ShowYourWork!]",
            default=user["affil"],
            style=style,
        ).ask()

        user["orcid"] = questionary.text(
            "ORCID:",
            qmark="[ShowYourWork!]",
            default=user["orcid"],
            style=style,
        ).ask()

        user["github"] = questionary.text(
            "GitHub Handle:",
            qmark="[ShowYourWork!]",
            default=user["github"],
            style=style,
        ).ask()

        # Save
        save_json(user, TEMP / "user" / "user.json")

        return user

    else:

        # Load from profile
        with open(TEMP / "user" / "user.json", "r") as f:
            user = json.load(f)

        return user
