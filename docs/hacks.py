import json
import os
import re
import subprocess
import sys
from pathlib import Path
from urllib.request import Request, urlopen

import jinja2
import requests
import yaml


def get_commit_count(repo, API_KEY):
    """
    Return the number of commits to a project.

    See https://stackoverflow.com/a/55873469

    """
    req = Request(f"https://api.github.com/repos/{repo}/commits?per_page=1")
    req.add_header("Accept", "application/vnd.github.v3+json")
    req.add_header("Authorization", f"token {API_KEY}")
    try:
        resp = urlopen(req)
        content = resp.read()
        content = json.loads(content)
        if len(content):
            commits = int(
                resp.headers["Link"]
                .split(",")[1]
                .split("per_page=1&page=")[1]
                .split(">")[0]
            )
        else:
            commits = 0
    except:
        commits = "N/A"
    return commits


def get_date(repo, API_KEY):
    """
    Get the timestamp when a repo was last pushed to.

    """
    req = Request(f"https://api.github.com/repos/{repo}")
    req.add_header("Accept", "application/vnd.github.v3+json")
    req.add_header("Authorization", f"token {API_KEY}")
    req.add_header("User-Agent", "request")
    try:
        content = urlopen(req).read()
        content = json.loads(content)
        date = content["pushed_at"]
    except:
        date = "??"
    return date


def get_version(repo, API_KEY):
    """
    Get the showyourwork SHA or version from the config file.

    """
    req = Request(f"https://api.github.com/repos/{repo}")
    req.add_header("Accept", "application/vnd.github.v3+json")
    req.add_header("Authorization", f"token {API_KEY}")
    try:
        resp = urlopen(req)
        content = resp.read()
        content = json.loads(content)
        branch = content.get("default_branch", "main")
        req = Request(
            f"https://raw.githubusercontent.com/{repo}/{branch}/showyourwork.yml"
        )
        content = urlopen(req).read()
        version = yaml.safe_load(content).get("version", "unknown")
    except:
        version = "unknown"
    return version


class RemoveSubmodulesHeading(StopIteration):
    pass


class RemoveContentsHeading(StopIteration):
    pass


# Snippets for each of the sections on the API page
SNAKEMAKE_SCRIPTS_DOCS = [
    """
These files are located in ``showyourwork/workflow/scripts`` and are executed
from the Snakemake rules defined in ``showyourwork/workflow/rules``.
They do the heavy lifting for the pipeline, including generating the article
graph in the preprocessing step, downloading and extracting Zenodo datasets,
and building the article PDF.
"""
]
SNAKEMAKE_RULES_DOCS = [
    """
These files are located in ``showyourwork/workflow/rules`` and contain the
main Snakemake rules for the pipeline. These rules either execute simple
shell commands (like running Python scripts to generate figures) or call
the scripts in ``showyourwork/workflow/scripts`` to perform more complex
pipeline tasks.
"""
]
SNAKEMAKE_SNAKEFILE_DOCS = [
    """
These are the two Snakefiles that control the pipeline. These Snakefiles
import all of the Snakemake rules defined in ``showyourwork/workflow/rules``,
parse the user config, ingest the user's custom Snakefile, etc.
"""
]
INTEGRATION_TEST_DOCS = [
    """
Integration tests are located in ``tests/integration`` at the root of the repo.
See below for information on how we test all the different moving parts in
``showyourwork``.
"""
]


# Add our module to the path
ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(1, str(ROOT / "showyourwork" / "workflow" / "scripts"))


# Generate the `projects.rst` page
API_KEY = os.getenv("GH_API_KEY", None)
if API_KEY is None:
    print("ERROR: Can't authenticate git. Unable to generate `projects.rst`.")
else:
    with open("projects.json", "r") as f:
        projects = json.load(f)
    for project in projects:
        projects[project]["date"] = projects[project].get("date", "")
        projects[project]["doi"] = projects[project].get("doi", "N/A")
        projects[project]["url"] = projects[project].get("url", "")
        projects[project]["version"] = get_version(project, API_KEY)
        projects[project]["commits"] = get_commit_count(project, API_KEY)
        projects[project]["date"] = get_date(project, API_KEY)
    fields = list(set([projects[project]["field"] for project in projects]))
    repos = sorted(projects.keys(), key=lambda item: projects[item]["date"])[
        ::-1
    ]
    count = {field: 0 for field in fields}
    for project in projects:
        count[projects[project]["field"]] += 1
    env = jinja2.Environment(loader=jinja2.FileSystemLoader("."))
    with open("projects.rst", "w") as f:
        print(
            env.get_template("projects.rst.jinja").render(
                projects=projects, fields=fields, repos=repos, count=count
            ),
            file=f,
        )


# Import the GitHub Action README
url = "https://raw.githubusercontent.com/showyourwork/showyourwork-action/main/README.rst"
try:
    r = requests.get(url, allow_redirects=True)
    content = r.content.decode()
    content = "The GitHub action\n=================\n\n" + "\n".join(
        content.split("\n")[10:]
    )
except:
    content = "The Github action\n=================\n\n"
    content += "Please visit `<https://github.com/showyourwork/showyourwork-action>`_."
    print("ERROR: Unable to generate `action.rst`.")
with open("action.rst", "w") as f:
    f.write(content)


# Create our API autodoc pages
for file in Path("api").rglob("*.rst"):
    file.unlink()
subprocess.run("sphinx-apidoc -feMT -d 5 -o api ../showyourwork", shell=True)
subprocess.run(
    "sphinx-apidoc -feM -d 5 -o api ../showyourwork/workflow/scripts",
    shell=True,
)


# Snakefile docs: ingest docstrings from the `.smk` files
rules = [
    str(file.name)
    for file in Path("../showyourwork/workflow/rules").glob("*.smk")
]
snakefiles = [
    str(file.name) for file in Path("../showyourwork/workflow").glob("*.smk")
]
for file in rules + snakefiles:
    if file in rules:
        with open(Path("../showyourwork/workflow/rules") / file, "r") as f:
            contents = f.read()
    else:
        with open(Path("../showyourwork/workflow") / file, "r") as f:
            contents = f.read()
    docstring = re.match('"""((?s).*?)"""', contents)
    if docstring:
        docstring = docstring.groups()[0]
    else:
        docstring = ""
    with open(f"api/{file}.rst", "w") as f:
        f.write(f"{file}\n{'=' * len(file)}\n\n{docstring}")

# Integration test docs
integration_tests = ["../integration_tests.rst"]


# Customize the API table of contents
with open("api/showyourwork.rst", "r") as f:
    lines = f.readlines()
lines = (
    [".. include:: ../api.rst\n"]
    + [
        "\n",
        ".. _api.module:",
        "\n",
        "\n",
        "The showyourwork module\n",
        "-----------------------\n",
    ]
    + lines[2:]
)
with open("api/modules.rst", "r") as f:
    new_lines = f.readlines()
Path("api/modules.rst").unlink()
lines += [
    "\n",
    ".. _api.scripts:",
    "\n",
    "\n",
    "Snakemake workflow scripts\n",
    "--------------------------\n",
    "\n",
    *SNAKEMAKE_SCRIPTS_DOCS,
    "\n",
] + new_lines[2:]
lines += [
    "\n",
    ".. _api.rules:",
    "\n",
    "\n",
    "Snakemake rules\n",
    "---------------\n\n",
    *SNAKEMAKE_RULES_DOCS,
    "\n",
    "\n",
    ".. toctree::\n",
    "   :maxdepth: 5\n",
    "\n   ",
    "\n   ".join(rules),
]
lines += [
    "\n\n",
    ".. _api.snakefiles:",
    "\n",
    "\n",
    "Snakefiles\n",
    "----------\n\n",
    *SNAKEMAKE_SNAKEFILE_DOCS,
    "\n",
    "\n",
    ".. toctree::\n",
    "   :maxdepth: 5\n",
    "\n   ",
    "\n   ".join(snakefiles),
]
lines += [
    "\n\n",
    ".. _api.integrationtests:",
    "\n",
    "\n",
    "Integration tests\n",
    "-----------------\n\n",
    *INTEGRATION_TEST_DOCS,
    "\n",
    "\n",
    ".. toctree::\n",
    "   :maxdepth: 5\n",
    "\n   ",
    "\n   ".join(integration_tests),
]
with open("api/showyourwork.rst", "w") as f:
    f.writelines(lines)


# Format each API rst file to get rid of unnecessary headings, etc.
for file in Path("api").rglob("*.rst"):
    with open(file, "r") as f:
        lines = f.readlines()
    lines[0] = (
        lines[0]
        .replace(" package", "")
        .replace(" module", "")
        .replace("showyourwork.", "")
    )
    while 1:
        try:
            subpackages = False
            for l, line in enumerate(lines):
                if line.startswith("Subpackages"):
                    subpackages = True
                    raise RemoveContentsHeading
                elif line.startswith("Submodules"):
                    if subpackages:
                        raise RemoveSubmodulesHeading
                    else:
                        raise RemoveContentsHeading
                elif line.startswith("Contents"):
                    raise RemoveContentsHeading
        except RemoveSubmodulesHeading:
            lines = lines[: l - 1] + lines[l + 6 :]
        except RemoveContentsHeading:
            lines = lines[:l] + lines[l + 3 :]
        else:
            break
    with open(file, "w") as f:
        f.writelines(lines)
