Developer API
=============

Welcome to the guts of |showyourwork|. The backend is fairly complicated given
that ``showyourwork`` is not just a Python package: it is simultanteously a
wrapper around Snakemake *and* a child process of Snakemake, as well as a command-line
tool, a cookiecutter template, a GitHub Action, and a conda environment manager. 
The code is wrapped up into a pip-installable Python package, but users should 
never have to actually import ``showyourwork`` into their scripts. Users should
only ever interact with the code via the command line interface, which spawns
child processes that themselves install ``showyourwork`` (since workflows *define*
specific versions of ``showyourwork`` to use).

But if you're reading this, you're probably interested in learning how ``showyourwork``
actually works, and possibly in tweaking it or adding new features. The documentation
pages below are divided into four sections. The first section contains documentation
for the Python package proper (:ref:`api.module`), followed by documentation 
for the scripts that do the heavy lifting during the article build step 
(:ref:`api.scripts`), the Snakemake rules that run those scripts (:ref:`api.rules`),
and the top-level Snakefiles that define the workflow (:ref:`api.snakefiles`).