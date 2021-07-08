import questionary
from pathlib import Path
import json


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
        "repo": {"url": "https://github.com/"},
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

    # Repository info
    defaults["repo"]["url"] = questionary.text(
        "URL:",
        qmark="[Repository]",
        default="https://github.com/{}/".format(defaults["author"]["github"]),
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
