The GitHub action
=================


The **showyourwork-action** runs on `GitHub Actions <https://github.com/features/actions>`_ to automatically build a `showyourwork <https://github.com/showyourwork/showyourwork>`_ article on the cloud every time changes are pushed to the remote repository. Under the hood, this action installs ``conda`` and ``showyourwork``, then runs the workflow to generate the article PDF, which it uploads to a separate branch on the remote. Importantly, everything is cached and all timestamps are preserved across builds, so this action will only re-run things that explicitly depend on the files that changed since the last time it ran.

This action is typically called in the workflow file ``.github/workflows/build.yml`` of a `showyourwork <https://github.com/showyourwork/showyourwork>`_ article repository. For more information on GitHub Actions workflow files, see `here <https://docs.github.com/en/actions/reference/workflow-syntax-for-github-actions>`_.

Inputs
------

The **showyourwork-action** accepts any of the following inputs, all of which are optional. These are provided using the ``with:`` directive in the ``showyourwork`` step of the ``.yml`` file, one per line (see the example below).

:code:`article-cache-number`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Optional** The **showyourwork-action** caches everything in your repository to speed up future builds. Sometimes, however, it's useful to clear the cache, such as when something breaks. This can be done by incrementing this number, which tells the action which version of the cache to load. Default: :code:`0`. Note that you can disable article caching by setting this variable to `null` or to an empty value.

:code:`build-tarball`
~~~~~~~~~~~~~~~~~~~~~

**Optional** Build a tarball for easy ArXiV submission? This tarball contains the article PDF, the rendered figures, and all the input files needed to compile the manuscript using a standard LaTeX compiler. The tarball is then pushed to the same branch as the article output (see ``force-push`` below). Default :code:`true`.

:code:`conda-cache-number`
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Optional** Bump this number to reset the :code:`conda` cache. The behavior is similar to that of ``article-cache-number`` above. Default: :code:`0`. Note that you can disable conda caching by setting this variable to `null` or to an empty value.

:code:`github-token`
~~~~~~~~~~~~~~~~~~~~

**Optional** A token for access to GitHub (e.g. :code:`secrets.GITHUB_TOKEN`). Do not set this value explicitly -- always use a secret! Default: :code:`${{ github.token }}` (usually set automatically).

:code:`output-branch-suffix`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Optional** Force-push output to branch :code:`<current-branch>-<output-branch-suffix>`? For example, if you've pushed a commit to the ``main`` branch, this action will by default compile your paper and force-push the output (the paper PDF as well as the ArXiV tarball, if enabled) to the branch ``main-pdf``. The *force* in *force-push* means this is not a typical ``git`` commit, as it will overwrite everything on that branch. This way, your repository won't get bloated over time with tons of committed output history. Default: :code:`pdf`.

Environment variables
---------------------

There are a few environment variables that may be needed on the ``showyourwork`` side. These include :code:`$ZENODO_TOKEN` and :code:`$SANDBOX_TOKEN` (API tokens that can be used to authenticate when publishing or downloading Zenodo/Zenodo Sandbox deposits) and :code:`$OVERLEAF_EMAIL` and :code:`$OVERLEAF_PASSWORD` (credentials for accessing and modifying an Overleaf project repository).
These should be provided through `Action secrets <https://docs.github.com/en/actions/security-guides/encrypted-secrets>`_ using the :code:`env:` directive (see the example below).

Concurrency
-----------

We recommend limiting the concurrency of **showyourwork-action** runs to one per branch. See `the docs <https://docs.github.com/en/actions/using-jobs/using-concurrency>`_ for details,
and check out the example below.

Example usage
-------------

Below is a complete example of a ``.github/workflows/build.yml`` file.

.. code-block:: yaml

  name: build

  on:
    push:
    workflow_dispatch:

  jobs:
    build:
      runs-on: ubuntu-latest
      name: Build the article PDF
      concurrency: showyourwork-${{ github.ref }}
      steps:

        - name: Checkout
          uses: actions/checkout@v2
          with:
            fetch-depth: 0

        - name: Build the article PDF
          id: build
          uses: showyourwork/showyourwork-action@v1
          env:
            ZENODO_TOKEN: ${{ secrets.ZENODO_TOKEN }}
            SANDBOX_TOKEN: ${{ secrets.SANDBOX_TOKEN }}
            OVERLEAF_EMAIL: ${{ secrets.OVERLEAF_EMAIL }}
            OVERLEAF_PASSWORD: ${{ secrets.OVERLEAF_PASSWORD }}