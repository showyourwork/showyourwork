.. _known_issues:

Known issues
============

This page collects issues that are known by the developers but they
are not yet fixed in the latest release.

Each entry reports also the original issue tracker from GitHub.

`Build on github is broken with 'pulp' has no attribute 'list_solvers' (#435) <https://github.com/showyourwork/showyourwork/issues/435>`_
-----------------------------------------------------------------------------------------------------------------------------------------

When using *showyourwork!*, along with the local version you installed on your computer
another one will be installed by GitHub for its server-side action
in order to compile the final document when you push your changes to the upstream repository.

If you see this error in your workflow logs you need to change what version of *showyourwork!*
is installed on the GitHub side.
This  can be directly done by editing the ``.github/workflows/build.yml`` file saved
in your local *showyourwork!* article directory.

For example, if you want GitHub to install and use the latest development version of *showyourwork!*,
then look for the following lines in the ``build.yml`` and ``build-pull-request.yml`` workflow files
and modify the ``Build the article PDF`` step like this,

.. code-block:: yaml

    - name: Build the article PDF
      id: build
      with:
        showyourwork-spec: git+https://github.com/showyourwork/showyourwork
      uses: showyourwork/showyourwork-action@v1
    