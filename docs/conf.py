# -- Project information -----------------------------------------------------

project = 'showyourwork'
copyright = '2021, Rodrigo Luger'
author = 'Rodrigo Luger'
release = '1.0.0'

# -- General configuration ---------------------------------------------------

extensions = []
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_book_theme'
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
