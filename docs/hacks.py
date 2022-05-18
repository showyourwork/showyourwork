import subprocess
import sys
from pathlib import Path
import re
import requests


# Add our module to the path
ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(1, str(ROOT / "showyourwork" / "workflow" / "scripts"))


class RemoveSubmodulesHeading(StopIteration):
    pass


class RemoveContentsHeading(StopIteration):
    pass


# Import the GitHub Action README
url = "https://raw.githubusercontent.com/showyourwork/showyourwork-action/main/README.rst"
r = requests.get(url, allow_redirects=True)
content = r.content.decode()
content = "The GitHub action\n=================\n\n" + "\n".join(
    content.split("\n")[10:]
)
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


# Customize the API table of contents
with open("api/showyourwork.rst", "r") as f:
    lines = f.readlines()
lines = (
    [".. include:: ../api.rst\n"]
    + [
        "\n",
        "The ``showyourwork`` module\n",
        "---------------------------\n",
    ]
    + lines[2:]
)
with open("api/modules.rst", "r") as f:
    new_lines = f.readlines()
Path("api/modules.rst").unlink()
lines += [
    "\n",
    "Snakemake workflow scripts\n",
    "--------------------------\n",
] + new_lines[2:]
lines += [
    "\n",
    "Snakemake rules\n",
    "---------------\n\n",
    ".. toctree::\n",
    "   :maxdepth: 5\n",
    "\n   ",
    "\n   ".join(rules),
]
lines += [
    "\n\n",
    "Snakefiles\n",
    "----------\n\n",
    ".. toctree::\n",
    "   :maxdepth: 5\n",
    "\n   ",
    "\n   ".join(snakefiles),
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