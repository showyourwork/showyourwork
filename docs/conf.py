import sys
from pathlib import Path
import ast
import builtins
import json
import jinja2
import os
import subprocess

# -- Hacks --------------------------------------------------------------------

# Add the workflow to the path
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "workflow"))
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "workflow" / "rules"))

# Add a global flag telling our workflow this is a sphinx run
builtins.__sphinx_docs_build__ = True

# Get list of projects that use showyourwork
if os.getenv("CI", "false") == "true":

    from get_repos import get_repos

    projects = get_repos()
    fields = list(
        set([projects[project]["field"] for project in projects])
        - set(["Uncategorized"])
    )
    repos = sorted(projects.keys(), key=lambda item: projects[item]["date"])[::-1]

    # Generate the `projects.rst` page
    env = jinja2.Environment(loader=jinja2.FileSystemLoader("."))
    with open("projects.rst", "w") as f:
        print(
            env.get_template("projects.rst.jinja").render(
                projects=projects, fields=fields, repos=repos
            ),
            file=f,
        )

# Get docstrings from scripts
scripts = Path(__file__).absolute().parents[1] / "workflow" / "scripts"
with open("scripts.rst", "w") as f:
    for script in scripts.glob("*.py"):
        with open(script, "r") as fin:
            tree = ast.parse(fin.read())
            docstring = ast.get_docstring(tree)
            print(script.name, file=f)
            print("^" * len(script.name), file=f)
            print(docstring, file=f)
            print("\n", file=f)

# Generate the `rules.rst` page
# NOTE: We don't set up Snakemake on RTD, so we need to run this
# locally and commit the `rules.rst` file periodically.
if not os.environ.get("READTHEDOCS") == "True":
    subprocess.check_call(
        ["snakemake", "-c1", "docstrings"], cwd=Path("..") / "workflow"
    )

# -- Project information -----------------------------------------------------

project = "showyourwork"
copyright = "2021, Rodrigo Luger"
author = "Rodrigo Luger"

# -- General configuration ---------------------------------------------------

extensions = ["sphinx.ext.autodoc"]
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
master_doc = "index"

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_book_theme"
html_copy_source = True
html_show_sourcelink = True
html_sourcelink_suffix = ""
html_title = "showyourwork"
html_logo = "_static/logo.png"
html_static_path = ["_static"]
html_css_files = []
html_theme_options = {
    "repository_url": "https://github.com/rodluger/showyourwork",
    "repository_branch": "main",
    "use_edit_page_button": True,
    "use_issues_button": True,
    "use_repository_button": True,
    "use_download_button": True,
    "logo_only": True,
    "use_fullscreen_button": False,
    "path_to_docs": "docs/",
}

# -- Extension settings ------------------------------------------------------

# autodoc
autoclass_content = "both"
autosummary_generate = True
autodoc_docstring_signature = True
autodoc_default_options = {"members": True}
