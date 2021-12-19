rule repo
^^^^^^^^^

Generates repository metadata (git url, branch name, commit sha).
Runs every time the article is generated.




rule aux_file
^^^^^^^^^^^^^

Copies auxiliary tex files to the ``src`` directory.




rule tectonic
^^^^^^^^^^^^^

Installs the latest version of tectonic. This rule isn't
used unless the user explicitly asks for ``tectonic_latest``
in the config file. We originally implemented this rule when
issues with ``archive.org`` caused the released versions of
``tectonic`` to fail, but an unreleased patch was functioning.
This rule is probably of very limited use in general.




rule class_name
^^^^^^^^^^^^^^^

Infers the name of the LaTeX document class for the article
by doing a regex on the TeX file. This helps ``showyourwork``
decide which class files to include in the build.




rule class_file
^^^^^^^^^^^^^^^

Copies the LaTeX class file to the ``src`` directory.




rule tmptexfile
^^^^^^^^^^^^^^^

Copies the manuscript to a temporary TeX file for the XML build step.




rule xml
^^^^^^^^

Generates the article XML tree. Specifically, builds the article
using ``tectonic``, but re-defines ``figure``, ``caption``, and ``label``
commands to print XML tags to a special log file. This way, we can
use LaTeX to construct a full XML tree of the document for us, without
any need for parsing the TeX file ourselves.
This XML tree is then used to determine relationships between the figure
scripts and the figure files.




rule script_info
^^^^^^^^^^^^^^^^

Builds the figure dependency tree from the XML file
generated in the ``xml`` rule. Stores it in a JSON
file in the temporary ``showyourwork`` directory.




rule metadata
^^^^^^^^^^^^^

Generates article metadata from the output of the ``repo``
and ``script_info`` rules. Saves it as a JSON in the temporary
``showyourwork`` directory.




rule figure
^^^^^^^^^^^

Generates a figure given a figure script and optional dependencies.
This is the workhorse of ``showyourwork``, which generates all the
figures in the article. Users can subclass this rule or override it
entirely to customize the build process.




rule texfile
^^^^^^^^^^^^

Copies the manuscript to a temporary TeX file for the final build.




rule stylesheet
^^^^^^^^^^^^^^^

Generates the ``showyourwork.sty`` LaTeX style sheet from a Jinja
template and copies it to the ``src`` directory.




rule dotgraph
^^^^^^^^^^^^^

Generate a .gv graph of the build process in the DOT language.




rule 14
^^^^^^^

Convert a figure file to a PNG thumbnail.




rule 15
^^^^^^^

Convert a figure file to a PNG thumbnail.




rule 16
^^^^^^^

Convert a figure file to a PNG thumbnail.




rule 17
^^^^^^^

Convert a figure file to a PNG thumbnail.




rule 18
^^^^^^^

Convert a figure file to a PNG thumbnail.




rule 19
^^^^^^^

Convert a figure file to a PNG thumbnail.




rule 20
^^^^^^^

Convert a figure file to a PNG thumbnail.




rule dag
^^^^^^^^

Generate a DAG (directed acyclic graph) of the build process.




rule pdf
^^^^^^^^

Builds the final article PDF. This is the final step in the article
build process.




rule arxiv
^^^^^^^^^^

Builds a tarball of the article PDF and all output for posting to the arXiv.




rule docstrings
^^^^^^^^^^^^^^^

Generates simple documentation for the rules in this workflow. This
rule is called when building the API documentation, but only when running
it locally (since we do not install ``Snakemake``, or ``conda`` for that
matter, on ``ReadTheDocs``). Therefore, we should build the documentation
locally from time to time and push changes to the ``rules.rst`` file
containing these docstrings so that our online docs are up to date.




