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


Generate an HTML build report:

.. code-block:: bash

    snakemake ms.pdf --report


Delete all output:

.. code-block:: bash

    snakemake -c1 ms.pdf --delete-all-output


Update ``showyourwork``:

.. code-block:: bash

    git submodule update --remote


Switch to specific version of ``showyourwork``:

    .. code-block:: bash

        pushd showyourwork
        git fetch --all --tags
        git checkout tags/v0.1.0
        popd
