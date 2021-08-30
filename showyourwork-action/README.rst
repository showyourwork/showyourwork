showyourwork-action
===================

The **showyourwork-action** makes it easy to build open-source, reproducible scientific articles from a LaTeX manuscript and Python figure scripts. This action is intended to be run on repositories generated using the `showyourwork <https://github.com/rodluger/showyourwork>`_ package.

This action generates all required figures from Python scripts in the :code:`figures/` directory of your repo and compiles the LaTeX manuscript :code:`tex/ms.tex` into the file :code:`ms.pdf`, which is uploaded as an artifact upon completion of the workflow. This file is also automatically force-pushed to a branch on your repository with the same name as the current branch plus the suffix :code:`-pdf`.

This action automatically caches all builds, so it won't re-run figure scripts that have not changed.

Inputs
------

:code:`action-path`
~~~~~~~~~~~~~~~~~~~

**Optional** Path to this action relative to the top level of the repo. Default: :code:`showyourwork/showyourwork-action`

:code:`article-cache-number`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Optional** Bump this number to reset the article cache. Default: :code:`0`

:code:`arxiv-tarball`
~~~~~~~~~~~~~~~~~~~~~

**Optional** Build a tarball for easy ArXiV submission? This will be pushed to the same branch as the article output. Default :code:`true`

:code:`arxiv-exclude`
~~~~~~~~~~~~~~~~~~~~~

**Optional** Files/paths to be excluded from the ArXiV tarball, one entry per line. Default :code:`**/*.py`, :code:`**/matplotlibrc`, :code:`**/.gitignore`

:code:`conda-cache-number`
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Optional** Bump this number to reset the :code:`conda` cache. Default: :code:`0`

:code:`conda-url`
~~~~~~~~~~~~~~~~~

**Optional** Exact url pointing to the :code:`conda` install script. Default: :code:`https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`

:code:`force-push`
~~~~~~~~~~~~~~~~~~

**Optional** Force-push output to branch :code:`<xxx>-pdf`, where :code:`<xxx>` is the current branch name? Default: :code:`true`

:code:`generate-report`
~~~~~~~~~~~~~~~~~~~~~~~

**Optional** Generate an HTML report of the build process and push it to GitHub Pages? Default: :code:`true`

:code:`gh-pages-branch`
~~~~~~~~~~~~~~~~~~~~~~~

**Optional** The branch serving GitHub Pages. Default: :code:`gh-pages`

:code:`github-token`
~~~~~~~~~~~~~~~~~~~~

**Optional** A token for access to GitHub (e.g. :code:`secrets.GITHUB_TOKEN`). Do not set this value explicitly -- always use a secret! Default: :code:`${{ github.token }}` (usually set automatically)

:code:`verbose`
~~~~~~~~~~~~~~~

**Optional** Enable verbose output and debug messages? Default: :code:`false`

Example usage
-------------

Below is a complete example that will automatically compile the article PDF on a repository with a :code:`showyourwork`-conforming layout. Note that it is recommended to check out the repo with `fetch-depth: 0` to enable output caching in between commits.

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
