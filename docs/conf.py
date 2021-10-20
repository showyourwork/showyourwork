import sys
from pathlib import Path
import builtins

# -- Hacks --------------------------------------------------------------------

# Add the workflow to the path
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "workflow"))
sys.path.insert(0, str(Path(__file__).absolute().parents[1] / "workflow" / "rules"))

# Add a global flag telling our workflow this is a sphinx run
builtins.__sphinx_docs_build__ = True


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
