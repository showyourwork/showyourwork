Contributor Guide
===================

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
users, so its helpful to create a separate environment for those. You should only ever
have to do this step once.

.. code-block:: bash

    conda env create -n sywdocs -f docs/environment.yml


During development
^^^^^^^^^^^^^^^^^^
The remaining steps will be run every time you work on the documentation.

1. Activate the docs environment

.. code-block:: text

    conda activate sywdocs


2. Go into the documentation folder

.. code-block:: text

    cd docs


3. Create the html from the reStructured Text (``.rst``) files

.. code-block:: text

    make html


4. Run a local file server of the docs website

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

Open Issues
---------------------------------
Another great way to contribute is to raise or open issues. Showyourwork issues
can be found `in the github repo <https://github.com/showyourwork/showyourwork/issues>`_. The
"contributions welcome" tag is a great place to start!
