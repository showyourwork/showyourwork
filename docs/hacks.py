import subprocess
import shutil
import sys
from pathlib import Path
import re
import requests


# Add our module to the path
ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(1, str(ROOT / "showyourwork" / "workflow" / "scripts"))


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
UNIT_TEST_DOCS = [
    """
Unit tests are located in ``tests`` at the root of the repo. See below for 
information on how we test all the different moving parts in ``showyourwork``.
"""
]


class RemoveSubmodulesHeading(StopIteration):
    pass


class RemoveContentsHeading(StopIteration):
    pass


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

# Unit test docs
unit_tests = ["../unit_tests.rst"]


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
    ".. _api.unittests:",
    "\n",
    "\n",
    "Unit tests\n",
    "----------\n\n",
    *UNIT_TEST_DOCS,
    "\n",
    "\n",
    ".. toctree::\n",
    "   :maxdepth: 5\n",
    "\n   ",
    "\n   ".join(unit_tests),
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