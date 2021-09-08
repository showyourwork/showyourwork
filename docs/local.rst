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
