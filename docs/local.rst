Local builds
============

Using ``make``
--------------

This page contains instructions on how to build your article locally; the same
applies to reproducing an article created by someone else using ``showyourwork``.
As within any GitHub repository, the first thing to do is clone it:

.. code-block:: bash

    git clone https://github.com/<user>/<repo>

Then, to build the paper, all you have to do is ``cd`` into the repository directory
and run

.. code-block:: bash

    make

This will initialize the ``showyourwork`` submodule and also install ``snakemake``
if you haven't already done these things manually (see below).

Note that if you cloned a third-party repository, and are looking to reproduce
their results, you may instead want to run

.. code-block:: bash

    make fast

which skips the generation step for any file that can instead be downloaded from
Zenodo. This can be useful if the paper depends on computationally expensive
steps whose output has been stored on Zenodo for quicker re-runs of the workflow.

You can also use the ``make`` command to produce individual figures, for example

.. code-block:: bash

    make src/figures/my_figure.pdf

or to delete all the output:

.. code-block:: bash

    make clean

Other commands include generating a directed acyclic graph (DAG) of the build process
(you'll need to have `graphviz <https://graphviz.org/download/>`_ installed):

.. code-block:: bash

    make dag

and generating an HTML build report for the workflow:

.. code-block:: bash

    make report

Under the hood, the ``Makefile`` calls ``Snakemake``, asking for one core by
default. You can change this behavior by providing a string of ``OPTIONS``, which
get passed to ``snakemake``:

.. code-block:: bash

    make ms.pdf OPTIONS="-c2"

for instance, to run the workflow on two cores. Type

.. code-block:: bash

    snakemake --help

for a list of all available options.


Using LaTeX Workshop in VSCode
------------------------------

If you edit and build your articles in VSCode, you will have to specify some
settings to enable VSCode to find the Makefile in your workspace root directory
rather then in the same directory as your manuscript file. Also, to view the PDF
in a parallel tab next to your manuscript file, you also have to tell LaTeX
Workshop where to find the PDF of the manuscript that ``showyourwork`` produces.
One solution that has worked for others is to create or edit a
workspace-specific settings file in ``.vscode/settings.json`` to add some custom
settings for LaTeX Workshop:

.. code-block:: json

    {
        "latex-workshop.latex.external.build.command": "make",
        "latex-workshop.latex.external.build.args": [],
        "latex-workshop.view.pdf.viewer": "tab",
        "latex-workshop.latex.outDir": "%WORKSPACE_FOLDER%"
    }

After this, you can use the ``LaTeX Workshop: Build LaTeX project`` command in
VSCode to build your manuscript file and have the PDF file auto-update in your
VSCode window.


Manual builds
-------------

While convenient, you don't need to use the ``Makefile`` to run
``showyourwork``. If you want to set up the repo manually, you should clone it
as follows

.. code-block:: bash

    git clone --recurse-submodules https://github.com/<user>/<repo>

or run

.. code-block:: bash

    git clone https://github.com/<user>/<repo>
    git submodule init
    git submodule update

to ensure the ``showyourwork`` submodule is downloaded and set up properly.

Next, if you don't already have them, install ``snakemake`` and ``jinja``:

.. code-block:: bash

    conda install -y -c defaults -c conda-forge -c bioconda mamba snakemake jinja

This step requires you to have the ``conda`` package manager
(click `here <https://www.anaconda.com/products/individual>`_ to download it).

Now, to build your paper, run

.. code-block:: bash

    snakemake -c1 --use-conda ms.pdf

from the top level of your repo.
This tells ``Snakemake`` to generate the file ``ms.pdf`` (your compiled article PDF)
on a single machine core (``-c1``) using the ``conda`` package manager.
The ``use-conda`` flag is imperative! But feel free to request more cores (``-c2``, ``-c3``, etc.)
if needed. You can also check out the `myriad options <https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_ you can pass to ``Snakemake``.

Some other useful commands:

- To generate a directed acyclic graph (DAG) of the build process, run

  .. code-block:: bash

      snakemake ms.pdf --dag | dot -Tpdf > dag.pdf


- To generate an HTML build report, run

  .. code-block:: bash

      snakemake ms.pdf --report


- To delete all output generated when running the ``ms.pdf`` rule, run

  .. code-block:: bash

      snakemake -c1 ms.pdf --delete-all-output
