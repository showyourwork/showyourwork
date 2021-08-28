Local builds
============

Install ``snakemake``:

.. code-block:: bash

    conda install -y -c defaults -c conda-forge -c bioconda mamba snakemake

Build the paper:

.. code-block:: bash

    snakemake -c1 --use-conda ms.pdf

Generate the directed acyclic graph (DAG):

.. code-block:: bash

    snakemake ms.pdf --dag | dot -Tpdf > dag.pdf

Delete all output:

.. code-block:: bash

    snakemake -c1 ms.pdf --delete-all-output
