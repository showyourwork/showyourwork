rule repo
^^^^^^^^^
Generate repository metadata (git url, branch name, commit sha).


rule aux_file
^^^^^^^^^^^^^
Copy auxiliary tex files to the ``src`` directory.


rule tectonic
^^^^^^^^^^^^^
Install the latest version of tectonic.


rule class_name
^^^^^^^^^^^^^^^
Infer the name of the LaTeX document class for the article.


rule class_file
^^^^^^^^^^^^^^^
Copy the LaTeX class file to the ``src`` directory.


rule tmptexfile
^^^^^^^^^^^^^^^
Copy the manuscript to a temporary TeX file for the XML build step.


rule xml
^^^^^^^^
Generate the article XML tree.


rule script_info
^^^^^^^^^^^^^^^^
Build the figure dependency tree.


rule metadata
^^^^^^^^^^^^^
Generate article metadata.


rule figure
^^^^^^^^^^^
Generate a figure given a figure script and optional dependencies.


rule texfile
^^^^^^^^^^^^
Copy the manuscript to a temporary TeX file for the final build.


rule stylesheet
^^^^^^^^^^^^^^^
Generating the ``showyourwork.sty`` LaTeX style sheet.


rule pdf
^^^^^^^^
Build the article PDF.


rule arxiv
^^^^^^^^^^
Build a tarball of the article PDF and all output for posting to the arXiv.


rule download_manual
^^^^^^^^^^^^^^^^^^^^
Download a figure dependency that was manually uploaded to Zenodo.


rule docstrings
^^^^^^^^^^^^^^^
Generate simple documentation for the rules in this workflow.


