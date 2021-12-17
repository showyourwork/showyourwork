from pathlib import Path

ROOT = Path(".").absolute().parents[1]
ROOT_URL = "https://github.com/rodluger/showyourwork-template/raw/main"

REQUIRED_DIRS = [
    ROOT / ".github",
    ROOT / ".github" / "workflows",
    ROOT / "src",
    ROOT / "src" / "figures",
]

REQUIRED_FILES = [
    ROOT / "README.md",
    ROOT / "LICENSE",
    ROOT / "Makefile",
    ROOT / "environment.yml",
    ROOT / "showyourwork.yml",
    ROOT / "Snakefile",
    ROOT / ".github" / "workflows" / "showyourwork.yml",
]

RECOMMENDED_DIRS = [ROOT / "src" / "static", ROOT / "src" / "data"]

RECOMMENDED_FILES = [
    ROOT / ".gitignore",
    ROOT / "src" / ".gitignore",
    ROOT / "src" / "figures" / ".gitignore",
    ROOT / "src" / "static" / ".gitignore",
    ROOT / "src" / "data" / ".gitignore",
    ROOT / "src" / "figures" / "matplotlibrc",
    ROOT / "CITATION.cff",
]


class Error:
    def __init__(self, *messages, level="error"):
        self.messages = messages
        self.level = level

    def __repr__(self):
        if self.level == "error":
            text = f"\033[1;31m[ ERROR ]\033[0m {self.messages[0]}"
        else:
            text = f"\033[1;32m[SUGGEST]\033[0m {self.messages[0]}"
        for message in self.messages[1:]:
            text += f"\n          {message}"
        return text


class lint:
    """
    Scans the user's repository for missing files or a wrong directory
    structure and recommends best-practice actions based on inspection
    of certain files.

    """

    def __init__(self):

        self.errors = []

        # Check the file structure
        self.check_tree()

        # Check the README.md
        self.check_readme()

        # Check the CITATION.cff
        self.check_citation()

        # Print errors
        if len(self.errors):
            for error in self.errors:
                print(error)
        else:
            print("\033[1;32mNo issues found. Keep up the great work!\033[0m")

    def check_tree(self):
        """
        Checks the repo tree for missing files / directories.

        """
        # Check for required folders
        for folder in REQUIRED_DIRS:
            if not folder.exists():
                self.errors.append(
                    Error(f"Missing required directory: {folder}/", level="error")
                )

        # Check for required files
        for file in REQUIRED_FILES:
            if not file.exists():
                self.errors.append(
                    Error(
                        f"Missing required file: {file}",
                        f"You may download it from: {ROOT_URL}/{file.relative_to(ROOT).as_posix()}",
                        level="error",
                    )
                )

        # Check for at least one TeX file
        if len(list((ROOT / "src").glob("*.tex"))) == 0:
            self.errors.append(Error(f"Missing a `*.tex` file in `src/`."))

        # Check for required folders
        for folder in RECOMMENDED_DIRS:
            if not folder.exists():
                self.errors.append(
                    Error(
                        f"Missing recommended directory: {folder}/", level="suggestion"
                    )
                )

        # Check for required files
        for file in RECOMMENDED_FILES:
            if not file.exists():
                self.errors.append(
                    Error(
                        f"Missing recommended file: {file}",
                        f"You may download it from: {ROOT_URL}/{file.relative_to(ROOT).as_posix()}",
                        level="suggestion",
                    )
                )

    def check_readme(self):
        """
        Checks the content of the README.

        """
        if not (ROOT / "README.md").exists():
            return

        with open(ROOT / "README.md", "r") as f:
            contents = f.read()

        # Check that there is a `make` instruction
        if "`make`" not in contents:
            if "```make```" not in contents.replace(" ", "").replace("\n", ""):
                self.errors.append(
                    Error(
                        "The README doesn't seem to have a `make` instruction.",
                        "We recommend including a section on how to reproduce your results!",
                        level="suggestion",
                    )
                )

        # Check that there is a `make fast` instruction
        if "`make fast`" not in contents:
            if "```make fast```" not in contents.replace(" ", "").replace("\n", ""):
                self.errors.append(
                    Error(
                        "The README doesn't seem to have a `make fast` instruction.",
                        "We recommend including a section on how to reproduce your results!",
                        level="suggestion",
                    )
                )

        # Check that the `showyourwork` logo has been removed
        if (
            '<img width = "450" src="https://raw.githubusercontent.com/rodluger/showyourwork/img/showyourwork.png" alt="showyourwork"/>'
            in contents
        ):
            self.errors.append(
                Error(
                    "It looks like you still have the `showyourwork` logo in your README.",
                    "While we are flattered, consider adding a header or logo showcasing your paper instead!",
                    level="suggestion",
                )
            )

    def check_citation(self):
        """
        Checks the content of the CITATION file.

        """
        if not (ROOT / "CITATION.cff").exists():
            return

        with open(ROOT / "CITATION.cff", "r") as f:
            contents = f.read()

        # Check that the user at least filled in their names
        if "Last Name" in contents:
            self.errors.append(
                Error(
                    "It looks like you haven't edited the CITATION.cff file.",
                    "Please take a moment to add the relevant info so others can cite your work.",
                    level="suggestion",
                )
            )
