Updating showyourwork
=====================

Update ``showyourwork``:

.. code-block:: bash

    git submodule update --remote


Switch to specific version of ``showyourwork``:

.. code-block:: bash

    pushd showyourwork
    git fetch --all --tags
    git checkout tags/v0.1.0
    popd
