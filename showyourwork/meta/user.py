from ..constants import *
from ..utils import save_json
import questionary
import subprocess
import json
import os


__all__ = ["get_user_metadata"]


def get_user_name():
    name = os.getenv("USER_EMAIL", "Author Name")
    if name == "Author Name":
        try:
            name = (
                subprocess.check_output(
                    ["git", "log", "--format='%an'", "HEAD^!"],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    cwd=USER,
                )
                .decode()
                .replace("\n", "")
            )
        except Exception as e:
            pass
    return name


def get_user_email():
    email = os.getenv("USER_EMAIL", "author@email")
    if email == "author@email":
        try:
            email = (
                subprocess.check_output(
                    ["git", "log", "--format='%ae'", "HEAD^!"],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    cwd=USER,
                )
                .decode()
                .replace("\n", "")
            )
        except Exception as e:
            pass
    return email


def get_user_github():
    return os.getenv("GITHUB_REPOSITORY", "unknown/").split("/")[0]


def get_user_metadata(clobber=True):

    if clobber or not (TEMP / "user" / "user.json").exists():

        if os.getenv("GITHUB_ACTIONS", None) is not None:

            user = {
                "name": "",
                "email": "",
                "affil": "",
                "orcid": "",
                "github": "",
            }

            # Use GitHub Actions settings
            user["name"] = get_user_name()
            user["email"] = get_user_email()
            user["affil"] = os.getenv("USER_AFFIL", "Author Affiliation")
            user["orcid"] = os.getenv("USER_ORCID", "0000-0000-0000-0000")
            user["github"] = get_user_github()

        else:

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
