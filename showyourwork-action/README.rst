Showyourwork GitHub action
==========================

The **showyourwork-action** runs on `GitHub Actions <https://github.com/features/actions>`_ to automatically build your article on the cloud every time you push changes to the repository. Under the hood, this action clones your repository, installs ``conda`` and ``Snakemake``, and runs your workflow to produce the figures and article PDF, which it then uploads to a separate branch on your repository. Importantly, everything is cached and all timestamps are preserved across builds, so this action will only re-run things that explicitly depend on the files that changed since the last time it ran.

This action is called in the workflow file ``.github/workflows/showyourwork.yml``, which is triggered on every pushed commit by default. Feel free to edit this file as needed, such as to include extra build steps, change the workflow triggers, or edit this action's settings (see below). For more information on GitHub Actions workflow files, see `here <https://docs.github.com/en/actions/reference/workflow-syntax-for-github-actions>`_.


Inputs
------

The **showyourwork-action** accepts any of the following inputs, all of which are optional. These are provided using the ``with:`` directive in the ``build`` step of the ``.yml`` file, one per line (see the example below).

:code:`action-path`
~~~~~~~~~~~~~~~~~~~

**Optional** Path to this action relative to the top level of the repo. You typically don't have to change this unless you've reorganized the layout of your repository (and you know what you're doing!). Default: :code:`showyourwork/showyourwork-action`.

:code:`article-cache-number`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Optional** The **showyourwork-action** caches everything in your repository to speed up future builds. Sometimes, however, it's useful to clear the cache, such as when something breaks. This can be done by incrementing this number, which tells the action which version of the cache to load. Default: :code:`0`.

:code:`arxiv-tarball`
~~~~~~~~~~~~~~~~~~~~~

**Optional** Build a tarball for easy ArXiV submission? This tarball contains the article PDF, the rendered figures, and all the input files needed to compile the manuscript using a standard LaTeX compiler. The tarball is then pushed to the same branch as the article output (see ``force-push`` below). Default :code:`true`.

:code:`arxiv-exclude`
~~~~~~~~~~~~~~~~~~~~~

**Optional** Files/paths to be excluded from the ArXiV tarball, one entry per line. If you have supplementary files or datasets that aren't explicitly required during the LaTeX compilation process, it's probably best to exclude them here! Default :code:`**/*.py`, :code:`**/matplotlibrc`, :code:`**/.gitignore`.

:code:`conda-cache-number`
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Optional** Bump this number to reset the :code:`conda` cache. The behavior is similar to that of ``article-cache-number`` above. Default: :code:`0`.

:code:`conda-url`
~~~~~~~~~~~~~~~~~

**Optional** Exact url pointing to the :code:`conda` install script. This always points to the latest ``conda`` installer, so you probably don't need to change this. Default: :code:`https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`.


:code:`github-token`
~~~~~~~~~~~~~~~~~~~~

**Optional** A token for access to GitHub (e.g. :code:`secrets.GITHUB_TOKEN`). Do not set this value explicitly -- always use a secret! Default: :code:`${{ github.token }}` (usually set automatically)

:code:`output-branch-suffix`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Optional** Force-push output to branch :code:`<current-branch>-<output-branch-suffix>`? For example, if you've pushed a commit to the ``main`` branch, this action will by default compile your paper and force-push the output (the paper PDF as well as the ArXiV tarball, if enabled) to the branch ``main-pdf``. The *force* in *force-push* means this is not a typical ``git`` commit, as it will overwrite everything on that branch. This way, your repository won't get bloated over time with tons of committed output history. Default: :code:`pdf`.

:code:`verbose`
~~~~~~~~~~~~~~~

**Optional** Enable verbose output and debug messages? Default: :code:`false`.


Example usage
-------------

Below is a complete example of a ``.github/workflows/showyourwork.yml`` file. Note that the repo should be checked out with ``fetch-depth: 0`` and ``submodules: recursive`` to ensure the action works as intended.

.. code-block:: yaml

    on: push
    jobs:
      showyourwork:
        runs-on: ubuntu-latest
        name: Build the article PDF
        steps:
          - name: Checkout
            uses: actions/checkout@v2
            with:
              fetch-depth: 0
              submodules: recursive

          - name: Build the article PDF
            id: build
            uses: ./showyourwork/showyourwork-action
            with:
              article-cache-number: 0
