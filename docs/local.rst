Local builds
============

Setup
-----

While ``showyourwork`` is mainly designed to build your article PDF on the cloud
via GitHub Actions, you'll likely want to build it locally as well. To get started
with that, clone your repository:

.. code-block:: bash

    git clone --recurse-submodules https://github.com/<user>/<repo>

or

.. code-block:: bash

    git clone https://github.com/<user>/<repo>
    git submodule init
    git submodule update

This makes sure the ``showyourwork`` submodule is downloaded and set up properly.

Next, if you don't already have it, install ``snakemake``:

.. code-block:: bash

    conda install -y -c defaults -c conda-forge -c bioconda mamba snakemake

This step requires you to have the ``conda`` package manager
(click `here <https://www.anaconda.com/products/individual>`_ to download it).

Build
-----

Now, to build your paper, navigate to the top level of your repository and run

.. code-block:: bash

    snakemake -c1 --use-conda ms.pdf

This tells ``Snakemake`` to generate the file ``ms.pdf`` (your compiled article PDF)
on a single machine core (``-c1``) using the ``conda`` package manager.
The ``use-conda`` flag is imperative! But feel free to request more cores (``-c2``, ``-c3``, etc.)
if needed. You can also check out the `myriad options <https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_ you can pass to ``Snakemake``.


Quick-and-dirty
---------------

As of version ``1.0.8``, the ``showyourwork-template`` repository comes with a ``Makefile``
that automates a lot of this stuff. After you clone your repo, you should be able to 
just run

.. code-block:: bash

    make ms.pdf

or, even simpler,

.. code-block:: bash

    make

to build your article. This will also install ``snakemake`` if you don't already
have it. You can also use this command to produce individual figures, for example

.. code-block:: bash

    make src/figures/my_figure.pdf

Under the hood, this is calling ``Snakemake`` with one core (``-c1``) and forcing
conda usage (``use-conda``). You can edit the ``Makefile`` variables to change these
options or add others.


Other
-----

To generate a directed acyclic graph (DAG) of the build process, run

.. code-block:: bash

    snakemake ms.pdf --dag | dot -Tpdf > dag.pdf


To generate an HTML build report, run

.. code-block:: bash

    snakemake ms.pdf --report


To delete all output generated when running the ``ms.pdf`` rule, run

.. code-block:: bash

    snakemake -c1 ms.pdf --delete-all-output

or simply

.. code-block:: bash

    make clean
