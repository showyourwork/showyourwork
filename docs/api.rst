Developer API
=============

Welcome to the guts of |showyourwork|. The backend is fairly complicated given
that ``showyourwork`` is not just a Python package: it is simultanteously a
wrapper around Snakemake *and* a child process of Snakemake, as well as a command-line
tool, a cookiecutter template, a GitHub Action, and a conda environment manager.
The code is wrapped up into a pip-installable Python package, but users should
never have to actually import ``showyourwork`` into their scripts. Users should
only ever interact with the code via the command line interface, which spawns
child processes that themselves import ``showyourwork``.

But if you're reading this, you're probably interested in learning how ``showyourwork``
actually works, and possibly in tweaking it or adding new features. The documentation
pages below are divided into four sections. The first section contains documentation
for the Python package proper (:ref:`api.module`), followed by documentation
for the scripts that do the heavy lifting during the article build step
(:ref:`api.scripts`), the Snakemake rules that run those scripts (:ref:`api.rules`),
and the top-level Snakefiles that define the workflow (:ref:`api.snakefiles`).

Contributing to the documentation
---------------------------------

A great place to start contributing is to documentation. In order to do this you'll
want to run the docs locally. Here are the steps for getting setup.

Initial setup
^^^^^^^^^^^^^

1. `Fork <https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_ the 
   |showyourwork| github repo and 
   `clone <https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository>`_ 
   your forked repo locally
2. Create a conda environment for the docs

Building the docs includes a few extra libraries that aren't required for |showyourwork|
users, so its helpful to create a seperate environment for those. You should only ever 
have to do this step once.

.. code-block:: bash

    cd docs
    conda env create -n sywdocs -f environment.yml


During development
^^^^^^^^^^^^^^^^^^
The remaining steps will be run every time you work on the documentation.

1. Activate the docs enviroment 

.. code-block:: text

    conda activate sywdocs


2. Create the html from the reStructured Text (``.rst``) files

.. code-block:: text

    make html


3. Run a local file server of the docs website

.. code-block:: text

    python3 -m http.server -d _build/html

From here you can go to a web browser and open up the docs using the
port shown in terminal. (Ex. ``localhost:8000``)

The last two steps will need to be run every time you make a change in the ``.rst``
files and want to view them. Because of this it can be nice to run them together using

.. code-block:: text

    make html && python3 -m http.server -d _build/html

That's the setup! Open documentation issues can be perused 
`here <https://github.com/showyourwork/showyourwork/issues?q=is%3Aissue+is%3Aopen+label%3A%22%3Amemo%3A+documentation%22>`_ .
