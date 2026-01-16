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

The simplest way to get an environment is to first create it

.. code-block:: bash

    conda create -n sywdocs python pip

activate it

.. code-block:: bash

    conda activate sywdocs

and install showyourwork for development

.. code-block:: bash

    python -m pip install -e ".[docs]"


Alternatively, you can create a conda environment

.. code-block:: bash

    conda create -n sywdocs -f docs/environment.yml

activate it

.. code-block:: bash

    conda activate sywdocs

and install showyourwork

.. code-block:: bash

    python -m pip install -e .

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

Contributing code
-----------------

1. start from an up-to-date branch of the main branch
   (name your branch something descriptive like 'fix-typo-in-docs' or 'add-new-feature')

.. code-block:: sh

    git switch main
    git pull upstream main # here 'upstream' is the original repo remote name, NOT your fork
    git switch -c my-feature-branch

2. make your changes trying to make small well-defined commits with descriptive messages
   (this will make it easier to review and understand the changes)
   and following the project's coding style

2. open a pull request against the main branch of the original repo

3. add a news fragment in the ``docs/changes`` directory

   This is required for all PRs. The news fragment should be a file named
   ``<pull_request_number>.<type>.rst``, ``<type>`` is one of the following:

   - ``bugfix``: for bug fixes
   - ``feature``: for new features
   - ``maintenance``: for maintenance tasks, like updating dependencies or documentation
   - ``api``: for changes to the public API
   - ``optimization``: for changes that improve performance or refactor code , but do not change user experience

Making a new release
--------------------

To generate the changelog before a new release,
execute the following command in the base directory of the project

.. code-block:: sh

    towncrier build --version=<VERSION NUMBER>

This will update and stage the CHANGES.rst file with the new release information,
removing the existing news fragments from the ``docs/changes`` directory.

You can then prepare a new PR with the title "Prepare release X.Y.Z" which should be merged
before proceeding the release.

TODO: Review instructions to actually make a new release.

Add new templates
-----------------

To add a new template to the list of available templates, you need to:

1. create a new template repository (please, make sure it works and it is up to date)
2. transfer the new template to the list of available templates under the _showyourwork_ organization
   (please, contact the administrators in case you do not have sufficient permissions)
3. add a link to the new template repository in the :ref:`templates`
   section of the documentation (start by copy-pasting an existing card entry and customize it for your use case)
