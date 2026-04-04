.. _install:

Installation
============

| In order to use |showyourwork|, you'll also need a working installation of the *conda* package manager.
| We recommend using `miniforge <https://conda-forge.org/download/>`_.

.. note::
    It is possible to run Snakemake without conda or mamba by
    passing the ``--no-conda`` flag (see :doc:`Command Line Interface <cli>`).
    However, we recommend using conda whenever possible to improve reproduciblity.
    When using showyourwork without conda, you must run it from an environment with
    all required software dependencies and a working
    `Tectonic installation <https://tectonic-typesetting.github.io>`_.

Create a new conda environment with the necessary dependencies:

.. code-block:: bash

    mamba create -n showyourwork pip
    mamba activate showyourwork

.. include:: dev_version_banner.rst

.. tab-set::

    .. tab-item:: Latest release

        .. code-block:: bash

            pip install showyourwork

    .. tab-item:: Development version

        .. code-block:: bash

            pip install git+https://github.com/showyourwork/showyourwork
