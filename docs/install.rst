.. _install:

Installation
============

| In order to use |showyourwork|, you'll also need a working installation of the *conda* package manager.
| We recommend using `miniforge <https://conda-forge.org/download/>`_.

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
