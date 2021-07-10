import questionary
from pathlib import Path
import json
import subprocess


# Paths
HOME = Path.home()
CORTEX = Path(__file__).parents[0].absolute()


def user_prompt():

    # Custom style
    style = questionary.Style([])

    # Load user defaults
    defaults = {
        "author": {
            "name": "",
            "email": "",
            "affil": "",
            "orcid": "",
            "github": "",
        },
        "repo": {"url": ""},
    }
    if (HOME / ".cortex").exists():
        with open(HOME / ".cortex", "r") as f:
            defaults.update(json.load(f))

    # Primary author info
    defaults["author"]["name"] = questionary.text(
        "Full Name:",
        qmark="[Author]",
        default=defaults["author"]["name"],
        style=style,
    ).ask()

    defaults["author"]["email"] = questionary.text(
        "Email:",
        qmark="[Author]",
        default=defaults["author"]["email"],
        style=style,
    ).ask()

    defaults["author"]["affil"] = questionary.text(
        "Affiliation:",
        qmark="[Author]",
        default=defaults["author"]["affil"],
        style=style,
    ).ask()

    defaults["author"]["orcid"] = questionary.text(
        "ORCID:",
        qmark="[Author]",
        default=defaults["author"]["orcid"],
        style=style,
    ).ask()

    defaults["author"]["github"] = questionary.text(
        "GitHub Handle:",
        qmark="[Author]",
        default=defaults["author"]["github"],
        style=style,
    ).ask()

    # Try to infer the GitHub url for this repo
    try:
        default_repo_url = subprocess.check_output(
            ["git", "config", "--get", "remote.origin.url"]
        ).decode()
        if default_repo_url.endswith("\n"):
            default_repo_url = default_repo_url[:-1]
        if default_repo_url.endswith(".git"):
            default_repo_url = default_repo_url[:-4]
    except Exception as e:
        default_repo_url = "https://github.com/{}".format(
            defaults["author"]["github"]
        )

    # Try to infer the current GitHub branch
    try:
        default_repo_branch = subprocess.check_output(
            ["git", "branch", "--show-current"]
        ).decode()
        if default_repo_branch.endswith("\n"):
            default_repo_branch = default_repo_branch[:-1]
    except Exception as e:
        default_repo_branch = "main"

    # Repository info
    url = questionary.text(
        "URL:",
        qmark="[Repository]",
        default=default_repo_url,
        style=style,
    ).ask()
    if url.endswith("/"):
        url = url[:-1]
    defaults["repo"]["url"] = url

    defaults["repo"]["branch"] = questionary.text(
        "Main branch:",
        qmark="[Repository]",
        default=default_repo_branch,
        style=style,
    ).ask()

    # Update user defaults
    with open(HOME / ".cortex", "w") as f:
        json.dump(defaults, f)

    # Update this repo's config
    with open(CORTEX / "data" / "meta.json", "w") as f:
        json.dump(defaults, f)


if __name__ == "__main__":
    user_prompt()
