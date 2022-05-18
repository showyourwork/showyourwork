Repository layout
=================

If you created your repository via `showyourwork setup <local>`_,
it should have the following overall layout:

.. raw:: html

      <style>
        /*
              https://codepen.io/asraven/pen/qbrQMX
        */
        .directory-list ul {
          margin-left: 10px;
          padding-left: 20px;
          border-left: 1px dashed #ddd;
        }
        .directory-list li {
          list-style: none;
          color: #888;
          font-size: 17px;
          font-style: normal;
          font-weight: normal;
        }
        .directory-list a {
          border-bottom: 1px solid transparent;
          color: rgba(var(--pst-color-link),1);
          text-decoration: none;
          transition: all 0.2s ease;
        }
        .directory-list a:hover {
          font-weight: bold;
        }
        .directory-list .folder > a {
          color: rgba(var(--pst-color-link),1);
        }
        .directory-list li:before {
          margin-right: 10px;
          content: "";
          height: 20px;
          vertical-align: middle;
          width: 20px;
          background-repeat: no-repeat;
          display: inline-block;
          /* file icon by default */
          background-image: url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><path fill='lightgrey' d='M85.714,42.857V87.5c0,1.487-0.521,2.752-1.562,3.794c-1.042,1.041-2.308,1.562-3.795,1.562H19.643 c-1.488,0-2.753-0.521-3.794-1.562c-1.042-1.042-1.562-2.307-1.562-3.794v-75c0-1.487,0.521-2.752,1.562-3.794 c1.041-1.041,2.306-1.562,3.794-1.562H50V37.5c0,1.488,0.521,2.753,1.562,3.795s2.307,1.562,3.795,1.562H85.714z M85.546,35.714 H57.143V7.311c3.05,0.558,5.505,1.767,7.366,3.627l17.41,17.411C83.78,30.209,84.989,32.665,85.546,35.714z' /></svg>");
          background-position: center 2px;
          background-size: 60% auto;
        }
        .directory-list li.folder:before {
          /* folder icon if folder class is specified */
          background-image: url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><path fill='lightblue' d='M96.429,37.5v39.286c0,3.423-1.228,6.361-3.684,8.817c-2.455,2.455-5.395,3.683-8.816,3.683H16.071 c-3.423,0-6.362-1.228-8.817-3.683c-2.456-2.456-3.683-5.395-3.683-8.817V23.214c0-3.422,1.228-6.362,3.683-8.817 c2.455-2.456,5.394-3.683,8.817-3.683h17.857c3.422,0,6.362,1.228,8.817,3.683c2.455,2.455,3.683,5.395,3.683,8.817V25h37.5 c3.422,0,6.361,1.228,8.816,3.683C95.201,31.138,96.429,34.078,96.429,37.5z' /></svg>");
          background-position: center top;
          background-size: 75% auto;
        }
      </style>

      <div class="box">
        <ul class="directory-list">
          <li class="folder">.github
            <ul>
              <li class="folder">workflow
                <ul>
                  <li><a href="#build">build.yml</a></li>
                  <li><a href="#build-pull-request">build-pull-request.yml</a></li>
                  <li><a href="#process-pull-request">process-pull-request.yml</a></li>
                </ul>
              </li>
            </ul>
          </li>
          <li class="folder"><a href="#src">src</a>
            <ul>
              <li class="folder"><a href="#data">data</a>
                <ul>
                  <li><a href="#gitignore">.gitignore</a></li>
                </ul>
              </li>
              <li class="folder"><a href="#scripts">scripts</a>
                <ul>
                  <li><a href="#gitignore">.gitignore</a></li>
                  <li><a href="#matplotlibrc">matplotlibrc</a></li>
                  <li><a href="#paths">paths.py</a></li>
                </ul>
              </li>
              <li class="folder"><a href="#static">static</a>
                <ul>
                  <li><a href="#gitignore">.gitignore</a></li>
                </ul>
              </li>
              <li class="folder"><a href="#tex">tex</a>
                <ul>
                  <li class="folder"><a href="#figures">figures</a>
                    <ul>
                      <li><a href="#gitignore">.gitignore</a></li>
                    </ul>
                  </li>
                  <li><a href="#gitignore">.gitignore</a></li>
                  <li><a href="#bibliography">bib.bib</a></li>
                  <li><a href="#manuscript">ms.tex</a></li>
                  <li><a href="#style">showyourwork.sty</a></li>
                </ul>
              </li>
            </ul>
          </li>
          <li><a href="#gitignore">.gitignore</a></li>
          <li><a href="#environment">environment.yml</a></li>
          <li><a href="#license">LICENSE</a></li>
          <li><a href="#readme">README.md</a></li>
          <li><a href="#config">showyourwork.yml</a></li>
          <li><a href="#snakefile">Snakefile</a></li>
          <li><a href="#zenodoconfig">zenodo.yml</a></li>
        </ul>
      </div>


Click on each of the files or folders in the directory tree above to read more about
them.

.. raw:: html

    <style>
      h2 {
        font-size: 21px !important;
      }
    </style>


.. _build:

The ``build.yml`` file
**********************

This is the configuration file for the workflow that builds your article on GitHub Actions.
It instructs GitHub Actions to build the article every time a commit is pushed to the
remote:

.. code-block:: yaml

    on:
      push:

The workflow tells GitHub to checkout the repository:

.. code-block:: yaml

    - name: Checkout
      uses: actions/checkout@v2
      with:
        fetch-depth: 0

and run the custom ``showyourwork-action`` to build the article:

.. code-block:: yaml

    - name: Build the article PDF
      id: build
      uses: showyourwork/showyourwork-action@v1

You can add other steps to this workflow or configure the action settings (see :doc:`action`), but most users shouldn't
have to tweak this file.
Check out the `GitHub documentation <https://docs.github.com/en/actions/reference/workflow-syntax-for-github-actions>`_ for
more information on configuring workflows.


.. _build-pull-request:

The ``build-pull-request.yml`` file
***********************************

The contents of this file are identical to those of :ref:`build`, except this workflow
is triggered only on pull requests. The reason for separating this out is that pull
request builds need a little bit of post-processing, which we can accomplish with
a special ``workflow_run`` trigger activated on runs of this workflow. See below
for details.


.. _process-pull-request:

The ``process-pull-request.yml`` file
*************************************

Pull request builds only get read access to the repository, so they can't upload the
article PDF anywhere. Instead, they generate a workflow artifact. This workflow
runs whenever a pull request build completes successfully. It downloads the build
artifact and force-pushes the article PDF to a special branch called (by default) 
``pull-request-<NUMBER>-pdf``, where ``<NUMBER>`` is the number of the pull request.
This workflow also posts a comment in the PR discussion with a link to the PDF
for convenience.


.. _src:

The ``src`` directory
*********************

This directory contains the source code for building your article. This includes
the LaTeX manuscript, the bibliography file, and the scripts
needed to produce all of the article figures. Figure scripts and auxiliary
code should be placed in the ``scripts`` directory; datasets should be 
programmatically generated or downloaded into the ``data`` directory; static
figures (if absolutely necessary!) should be placed in the ``static``
directory; and TeX files should be placed in the ``tex`` directory. 
See below for details.


.. _data:

The ``data`` directory
************************

This directory is included in the template as a convenience. It is meant to
house temporary (non-tracked) datasets, such as those downloaded from Zenodo
or programmatically generated by a pipeline script.
By default, nothing in this directory is tracked by ``git``.


.. _gitignore:

The ``.gitignore`` files
************************

The ``.gitignore`` files prevent you from committing certain kinds of files.
You can add entries to these files, but you shouldn't have to remove any.
In general, you should never commit figures (``.pdf``, ``.png``, ``.tiff``, etc),
LaTeX temporaries (``.log``, ``.aux``, ``.blg``, etc), or any kind of output.
In some cases, it might make sense to include one of these files (say, a ``.png``
photograph that can't be generated programatically from a script). To commit
something that's disallowed by a ``.gitignore`` file, just use the ``-f`` or ``--force``
option when adding the file to ``git``.


.. _scripts:

The ``scripts`` directory
*************************

This directory should contain all of the scripts needed to produce the figures
and/or datasets used in the generation of
your article.


.. _matplotlibrc:

The ``matplotlibrc`` config
***************************

This is a ``matplotlib`` configuration file
(see `the matplotlib docs <https://matplotlib.org/stable/tutorials/introductory/customizing.html>`_).
By default, it contains only a single instruction:

.. code-block::

    backend: agg

which tells ``matplotlib`` to render figures using the non-interactive backend
``agg``. You can change this and add additional config options, but note that
if the figures are generated using an interactive backend, this will likely
cause the GitHub Action to fail!


.. _paths:

The ``paths.py`` file
*********************

This file is included as a convenience to make it easy for scripts to load and save files
from/to the correct directories. It defines variables such as ```data``, ``figures``,
``scripts``, etc.; these are ``pathlib.Path`` instances corresponding to the absolute
path to the directories of the same name in your repository.
As an example, the following code

.. code-block:: python

    import paths
    import numpy as np
    import matplotlib.pyplot as plt

    data = np.loadtxt(paths.data / "dataset.txt")
    fig = plt.figure()
    plt.plot(data)
    fig.savefig(paths.figures / "figure.pdf")
    
loads a dataset called ``dataset.txt`` from the :ref:`data<data>` directory, plots it,
and saves the figure as ``figure.pdf`` in the :ref:`figures<figures>` directory.
All paths declared in the ``paths`` module are absolute, so the above script will
work regardless of what directory it is executed from.


.. _static:

The ``static`` directory
************************

This directory is meant to house figure files that can't be generated from
scripts, such as photos, flowcharts, reproductions of figures in other papers, etc.
If you place your figure in here, |showyourwork| will know not to try to
generate it from any script.


.. _tex:

The ``tex`` directory
*********************

This is the directory containing the TeX manuscript, your bibliography, and any
other auxiliary files used in building the article PDF, such as custom class
files or TeX style sheets. The contents of this folder can be synced to/from
an Overleaf project (see :doc:`overleaf`). The subfolder ``figures`` (see below)
contains the actual figure files to be included in the article build.


.. _figures:

The ``figures`` directory
*************************

This directory contains the figure output (the files generated by the figure
scripts) that gets included in your final article PDF. The ``.gitignore`` file
in this directory prevents anything in it from being tracked by ``git``, as
all figures should be programmatically generated on the fly.


.. _bibliography:

The ``bib.bib`` bibliography
****************************

This is an optional LaTeX/BibTeX bibliography file. Feel free to delete or rename
it if needed.


.. _manuscript:

The ``ms.tex`` manuscript
*************************

This is your manuscript file. By default, it's a `AASTeX v6.3.1 <https://journals.aas.org/aastexguide/>`_
article file with a placeholder title, abstract, and introduction. Feel free to
change the article class to whatever you'd like, although you may have to include
(and commit) the ``.sty`` stylesheet if it's not in the standard TeXLive distribution.
You can also import whatever packages you want in ``ms.tex`` -- the ``tectonic``
typesetting system will automatically download and install them when building
your article. Note that you can also rename this file to something else, as
long as you edit the corresponding setting (see :doc:`config`).


.. _style:

The ``showyourwork.sty`` file
*****************************

This is the |showyourwork| TeX style sheet, which you should always include in
your manuscript:

.. code-block:: tex

    \usepackage{showyourwork}

If you peek inside, there's not much there: it's a placeholder stylesheet that
imports all the |showyourwork| functionality from elsewhere if the article
is built as part of the |showyourwork| workflow. If you compile your article
with a standard TeX compiler (such as ``pdflatex`` or ``tectonic``), things will
still work, but you won't benefit from any of the showyourwork functionality.

.. _environment:

The ``environment.yml`` file
****************************

The ``environment.yml`` file specifies all of the ``conda`` packages needed to
build your article. You can read more about environment files
`in the conda docs <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.
By default, only the bare minimum specs are included
(e.g., ``numpy`` and ``matplotlib``). Feel free to manually add to this list
(noting that packages that can only be installed via ``pip`` should be placed in
the ``pip`` section). It's recommended to either pin a specific version
(i.e., ``matplotlib==3.3.4``) or specify a minimum version
(i.e., ``matplotlib>=3.0.0``) for your packages. Just be aware that overconstrained
requirements may break on other platforms
(see `this post <https://stackoverflow.com/questions/39280638/how-to-share-conda-environments-across-platforms>`_),
so you should probably only pin the direct dependencies of your project.
If you alread have a ``conda`` environment for your project, you can export
these direct dependencies -- the ones that you explicitly installed in the enviornment --
by running

.. code-block:: bash

    conda env export --from-history | grep -v "^prefix: " > environment.yml

The ``grep`` command removes the line in the environment file with the absolute path
to your ``conda`` environment, which probably won't be useful to anyone else running
your code!


.. _license:

The ``LICENSE`` file
********************

The ``LICENSE`` included in your repository is by default the MIT open-source
license. Feel free to change this to whatever license you prefer, although we
strongly recommend you to keep everything open source and free for others to
modify and adapt into their own work!


.. _readme:

The ``README.md`` file
**********************

You should include a description of your repository here. Keep the badges at
the top, as these provide easy access to the compiled article and the build logs.
Feel free to remove or change the logo, though!


.. _config:

The ``showyourwork.yml`` config file
************************************

This is the `Snakemake config file <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html>`_
for |showyourwork|, where you can customize several aspects of the build. For
detailed information on this file, see :doc:`the showyourwork.yml file <config>`.


.. _snakefile:

The ``Snakefile`` workflow
**************************

The ``Snakefile`` contains custom instructions to build your article
from the files in your repository. If you're not familiar with the Snakemake
workflow management system, read up on it `here <https://snakemake.readthedocs.io>`_.
By default, the ``Snakefile`` is empty: |showyourwork| takes care of everything
for you. For custom workflows, you can add rules to your ``Snakefile``, such
as instructions on how to build custom figures, to download datasets, etc; see
:doc:`custom`.
Note, finally, that this file is written in a language that's a straightforward
superset of Python, so any regular Python commands and syntax is valid in it.


.. _zenodoconfig:

The ``zenodo.yml`` config file
******************************

This is a config file used to store Zenodo caching information. This file is
entirely managed by |showyourwork|, so users should not edit it manually.