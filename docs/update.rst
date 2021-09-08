Updating showyourwork
=====================

If you would like to update to the latest version of ``showyourwork``,
check out the :doc:`Changelog <changelog>` to take a look at what's changed.
In most cases all you need to do is to update the ``showyourwork`` submodule:

.. code-block:: bash

    git submodule update --remote

Then add, commit, and push your changes to GitHub. Note that future versions of
``showyourwork`` may change the repository layout or the workflow syntax, so
be sure to read the docs carefully when upgrading!

You can also switch to a specific version of ``showyourwork`` (such as to roll it
back in case the upgrade broke something!). First, check the
`releases page on GitHub <https://github.com/rodluger/showyourwork/releases>`_
to see what versions are available, then simply checkout the desired one:

.. code-block:: bash

    pushd showyourwork
    git fetch --all --tags
    git checkout tags/v0.1.0
    popd
