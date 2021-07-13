from cortex_version import __version__
import glob
from pathlib import Path
import shutil
import subprocess
import jinja2
import json
import xml.etree.ElementTree as ET
import numpy as np
import warnings
import questionary
import os


# Path
ROOT = Path(__file__).parents[1].absolute()


# Git error codes
ctxScriptDoesNotExist = 3
ctxScriptNotVersionControlled = 2
ctxScriptHasUncommittedChanges = 1
ctxScriptUpToDate = 0
ERROR_CODES = {
    ctxScriptDoesNotExist: "ctxScriptDoesNotExist",
    ctxScriptNotVersionControlled: "ctxScriptNotVersionControlled",
    ctxScriptHasUncommittedChanges: "ctxScriptHasUncommittedChanges",
    ctxScriptUpToDate: "ctxScriptUpToDate",
}


# Special Jinja env for LaTeX
JINJA_ENV = jinja2.Environment(
    block_start_string="((*",
    block_end_string="*))",
    variable_start_string="((-",
    variable_end_string="-))",
    comment_start_string="((=",
    comment_end_string="=))",
    trim_blocks=True,
    autoescape=False,
    loader=jinja2.FileSystemLoader(ROOT / ".cortex" / "templates"),
)


class TemporaryTexFiles:
    def __init__(self, **kwargs):
        with open(ROOT / ".cortex" / "tex" / "cortex.sty", "w") as f:
            cortex_sty = JINJA_ENV.get_template("cortex.sty.jinja").render(
                **kwargs
            )
            print(cortex_sty, file=f)

    def __enter__(self):
        self.files = glob.glob(str(ROOT / ".cortex" / "tex" / "*"))
        for file in self.files:
            shutil.copy(file, ROOT / "tex")

    def __exit__(self, exc_type, exc_value, exc_traceback):
        files = [str(ROOT / "tex" / Path(file).name) for file in self.files]
        for suff in ["*.log", "*.blg", "*.xml"]:
            for f in glob.glob(str(ROOT / "tex" / suff)):
                files.append(f)
        for file in files:
            os.remove(file)


def get_git_status(path, script):
    # Check if the file exists
    if not (ROOT / path / script).exists():
        error_code = ctxScriptDoesNotExist
    else:
        # Check if the file is version controlled
        try:
            status = (
                subprocess.check_output(
                    ["git", "ls-files", "--error-unmatch", script],
                    cwd=ROOT / path,
                    stderr=subprocess.DEVNULL,
                )
                .decode()
                .replace("\n", "")
            )
        except Exception as e:
            # File is not version controlled
            error_code = ctxScriptNotVersionControlled
        else:
            # Check if the file has uncommitted changes
            try:
                status = (
                    subprocess.check_output(
                        ["git", "status", "-s", script],
                        cwd=ROOT / path,
                        stderr=subprocess.DEVNULL,
                    )
                    .decode()
                    .replace("\n", "")
                )
                if len(status):
                    raise Exception("Uncommitted changes!")
            except Exception as e:
                # Assume uncommited changes
                error_code = ctxScriptHasUncommittedChanges
            else:
                # File is good!
                error_code = ctxScriptUpToDate
    return error_code


def get_equation_labels(equation):
    return [label.text for label in equation.findall("LABEL")]


def get_equation_script(label):
    return "test_{}.py".format(label)


def check_figure_format(figure):

    # Get all figure elements
    elements = list(figure)
    captions = figure.findall("CAPTION")
    labels = figure.findall("LABEL")

    # Check that figure labels aren't nested inside captions
    for caption in captions:
        caption_labels = caption.findall("LABEL")
        if len(caption_labels):
            raise ValueError(
                "Label `{}` should not be nested within the figure caption".format(
                    caption_labels[0].text
                )
            )

    # The label must always come after the figure caption
    if len(captions):

        # Index of last caption
        caption_idx = (
            len(elements)
            - 1
            - np.argmax(
                [element.tag == "CAPTION" for element in elements[::-1]]
            )
        )

        if len(labels):

            # Index of first label
            label_idx = np.argmax(
                [element.tag == "LABEL" for element in elements]
            )

            if label_idx < caption_idx:
                raise ValueError(
                    "Figure label `{}` must come after the caption.".format(
                        labels[0].text
                    )
                )

    # Check that there is exactly one label
    if len(labels) >= 2:
        raise ValueError(
            "A figure has multiple labels: `{}`".format(", ".join(labels))
        )
    elif len(labels) == 0:
        # Unable to infer script
        warnings.warn("There is a figure without a label.")


def get_figure_label(figure):
    return figure.findall("LABEL")[0].text


def get_figure_script(label):
    return "{}.py".format(label)


def get_figure_files(figure):
    return [graphic.text for graphic in figure.findall("GRAPHICS")]


def save_json(contents, file):
    """
    Saves a dictionary in JSON format to disk only if the contents
    of the file would be changed (or if the file doesn't exist).

    """
    current = json.dumps(contents, indent=4)
    if Path(file).exists():
        with open(file, "r") as f:
            existing = json.load(f)
        if existing != current:
            with open(file, "w") as f:
                print(current, file=f)
    else:
        with open(file, "w") as f:
            print(current, file=f)


def get_user_metadata(clobber=False):

    if clobber or not (ROOT / ".cortex" / "data" / "user.json").exists():

        # Custom style
        style = questionary.Style([])

        # Load user defaults
        user = {
            "author": {
                "name": "",
                "email": "",
                "affil": "",
                "orcid": "",
                "github": "",
            },
            "repo": {"url": ""},
        }
        if (Path.home() / ".cortex").exists():
            with open(Path.home() / ".cortex", "r") as f:
                user.update(json.load(f))

        # Primary author info
        user["author"]["name"] = questionary.text(
            "Full Name:",
            qmark="[Author]",
            default=user["author"]["name"],
            style=style,
        ).ask()

        user["author"]["email"] = questionary.text(
            "Email:",
            qmark="[Author]",
            default=user["author"]["email"],
            style=style,
        ).ask()

        user["author"]["affil"] = questionary.text(
            "Affiliation:",
            qmark="[Author]",
            default=user["author"]["affil"],
            style=style,
        ).ask()

        user["author"]["orcid"] = questionary.text(
            "ORCID:",
            qmark="[Author]",
            default=user["author"]["orcid"],
            style=style,
        ).ask()

        user["author"]["github"] = questionary.text(
            "GitHub Handle:",
            qmark="[Author]",
            default=user["author"]["github"],
            style=style,
        ).ask()

        # Try to infer the GitHub url for this repo
        try:
            default_repo_url = subprocess.check_output(
                ["git", "config", "--get", "remote.origin.url"],
                stderr=subprocess.DEVNULL,
            ).decode()
            if default_repo_url.endswith("\n"):
                default_repo_url = default_repo_url[:-1]
            if default_repo_url.endswith(".git"):
                default_repo_url = default_repo_url[:-4]
        except Exception as e:
            default_repo_url = "https://github.com/{}".format(
                user["author"]["github"]
            )

        # Try to infer the current GitHub branch
        try:
            default_repo_branch = subprocess.check_output(
                ["git", "branch", "--show-current"], stderr=subprocess.DEVNULL
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
        user["repo"]["url"] = url

        user["repo"]["branch"] = questionary.text(
            "Main branch:",
            qmark="[Repository]",
            default=default_repo_branch,
            style=style,
        ).ask()

        # Update user defaults
        save_json(user, Path.home() / ".cortex")

        # Update this repo's config
        save_json(user, ROOT / ".cortex" / "data" / "user.json")

        return user

    else:

        with open(ROOT / ".cortex" / "data" / "user.json", "r") as f:
            defaults = json.load(f)

        return defaults


def get_script_metadata(clobber=False):

    if clobber or not (ROOT / ".cortex" / "data" / "user.json").exists():

        # Generate the XML file
        with TemporaryTexFiles(gen_tree=True):
            subprocess.check_output(
                ["tectonic", "--keep-logs", "-r", "0", "ms.tex"],
                cwd=str(ROOT / "tex"),
            )
            os.remove(ROOT / "tex" / "ms.pdf")
            root = ET.parse(ROOT / "tex" / "cortex.xml").getroot()

        # Parse figures
        figures = {}
        for figure in root.findall("FIGURE"):
            check_figure_format(figure)
            label = get_figure_label(figure)
            script = get_figure_script(label)
            files = get_figure_files(figure)
            status = get_git_status("figures", script)
            figures[label] = {
                "script": script,
                "files": files,
                "status": status,
            }

        # Parse equations
        equations = {}
        for equation in root.findall("ALIGN"):
            labels = get_equation_labels(equation)
            for label in labels:
                script = get_equation_script(label)
                status = get_git_status("tests", script)
                equations[label] = {"script": script, "status": status}

        # Store as JSON
        scripts = {"figures": figures, "equations": equations}
        save_json(scripts, ROOT / ".cortex" / "data" / "scripts.json")

        return scripts

    else:

        with open(ROOT / ".cortex" / "data" / "scripts.json", "r") as f:
            scripts = json.load(f)

        return scripts


def get_metadata(clobber=False):

    if clobber or not (ROOT / ".cortex" / "data" / "meta.json").exists():

        # Load the metadata (no clobber)
        user = get_user_metadata()
        scripts = get_script_metadata()

        # Get the current git hash
        meta = dict(user)
        try:
            meta["sha"] = (
                subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
                )
                .decode()
                .replace("\n", "")
            )
        except:
            meta["sha"] = meta["repo"]["branch"]

        # Miscellaneous
        meta["version"] = __version__
        meta["gen_tree"] = False

        # Figure and equation metadata
        status = 0
        meta["labels"] = {}
        for label, value in scripts["figures"].items():
            meta["labels"]["{}_script".format(label)] = value["script"]
            meta["labels"]["{}_status".format(label)] = ERROR_CODES[
                value["status"]
            ]
            status = max(status, value["status"])
        for label, value in scripts["equations"].items():
            meta["labels"]["{}_script".format(label)] = value["script"]
            meta["labels"]["{}_status".format(label)] = ERROR_CODES[
                value["status"]
            ]
            status = max(status, value["status"])
        meta["status"] = ERROR_CODES[status]

        # Store as JSON
        save_json(meta, ROOT / ".cortex" / "data" / "meta.json")

        return meta

    else:

        with open(ROOT / ".cortex" / "data" / "meta.json", "r") as f:
            meta = json.load(f)

        return meta


def gen_pdf():
    with TemporaryTexFiles(**get_metadata()):
        subprocess.check_output(
            ["tectonic", "--keep-logs", "ms.tex"],
            cwd=str(ROOT / "tex"),
        )
        shutil.move(ROOT / "tex" / "ms.pdf", ROOT)
