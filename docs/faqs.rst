FAQs
====

.. note:: This page is still under construction. More FAQs will be added shortly!

Permissions errors in GitHub Actions
------------------------------------

If you try to build a |showyourwork| article from a repository fork
on GitHub Actions, you may run into the following error when the action
attempts to push the results to the ``-pdf`` branch:

.. code-block:: text

    Uploading output
    /tmp/tmp.KORuhtnUA7
    Switched to a new branch 'main-pdf'
    [main-pdf (root-commit) 224ecfd] force-push article output
    2 files changed, 0 insertions(+), 0 deletions(-)
    create mode 100644 arxiv.tar.gz
    create mode 100644
    remote: Permission to $USER/$REPO.git denied to github-actions[bot].
    fatal: unable to access 'https://github.com/$USER/$REPO/': The requested URL returned error: 403

This happens because the default GitHub Actions permissions for the ``GITHUB_TOKEN``
secret are set to ``restricted`` for repository forks. To get the build to work,
go to

.. raw:: html

    <pre>
    https://github.com/<span class="text-highlight">$USER/$REPO</span>/settings/actions
    </pre>

and change the permissions to ``permissive``:

.. image:: _static/workflow_permissions.png
   :width: 60%
   :align: center


Rendering LaTeX in matplotlib
-----------------------------

When plotting with ``matplotlib``, if you run into errors that look like

.. code-block:: text

    FileNotFoundError: [Errno 2] No such file or directory: 'latex'

or

.. code-block:: text

    RuntimeError: Failed to process string with tex because latex could not be found

you are probably missing a proper ``latex`` installation. Recall that |showyourwork|
uses ``tectonic`` to build your article, which is not compatible with ``matplotlib``.
Instead, you'll have to install a separate TeX distribution, such as TeXLive or MiKTeX.
The same applies to runs on GitHub actions.

The simplest workaround is to disable LaTeX rendering in ``matplotlib``:

.. code-block:: python

    import matplotlib.pyplot as plt
    plt.rcParams.update({"text.usetex": False})

Math-mode strings can still be parsed using the built-in ``matplotlib`` renderer,
and in most cases this will do what you need. In some cases, however, the built-in
renderer may not cut it. If you really need a proper LaTeX installation, you'll
have to do a bit of extra work to get your build to pass on GitHub Actions.
First, add the following step to the ``build.yml`` and ``build-pull-request.yml`` 
workflows in your ``.github/workflows`` folder, just before the |showyourwork| 
``build`` step:

.. code-block:: yaml

    - name: Install TinyTex for matplotlib LaTeX rendering
      id: tinytex
      shell: bash -l {0}
      run: |
        wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh" | sh 
        sudo ~/bin/tlmgr install type1cm cm-super

This will install `TinyTex <https://yihui.org/tinytex/>`_, a
very lightweight TeX distribution that should provide everything you need. Note
that this step also installs the ``type1cm`` and ``cm-super`` LaTeX packages,
which may be required by ``matplotlib``. You can specify additional packages
in the same line if needed.

Then, in order for ``matplotlib`` to execute ``latex``, the ``~/bin`` path needs to 
be in the system ``$PATH``. This variable gets overwritten when running scripts inside isolated
``conda`` environments (as |showyourwork| does), so you'll need to add ``~\bin``
to the ``$PATH`` *within* your Python script. Therefore, add the
following bit of boilerplate to the top of any scripts that require LaTeX parsing:

.. code-block:: python

    import os
    from pathlib import Path
    os.environ["PATH"] += os.pathsep + str(Path.home() / "bin")

To save some typing, you could instead add this boilerplate to the 
``src/scripts/paths.py`` file so that
these commands get executed whenever that file is imported into your scripts.


Using LaTeX Workshop in VSCode
------------------------------

If you edit and build your articles in VSCode, you need to specify some settings so that VSCode knows to use |showyourwork| to build your document.
You can do this by creating (or editing) a workspace-specific settings file, ``.vscode/settings.json``, in the root directory of your repo.
At minimum, you should add the following lines:

.. code-block:: json

    {
        "latex-workshop.latex.external.build.command": "showyourwork",
        "latex-workshop.latex.external.build.args": [],
        "latex-workshop.latex.outDir": "%WORKSPACE_FOLDER%",
        "latex-workshop.view.pdf.viewer": "tab",
    }

This enables you to build the document using ``LaTeX Workshop: Build LaTeX project`` in the command palette.
Note that the final line tells LaTeX Workshop to open your article pdf in a VSCode tab.
Feel free to change ``tab`` to ``browser`` if you would rather LaTeX Workshop open your article in a browser tab.


If you also want to use LaTeX Workshop's AutoBuild on save (or on file change), you can add the following lines to the settings file

.. code-block:: json

    {
        "latex-workshop.latex.recipe.default": "showyourwork",
        "latex-workshop.latex.recipes": [
            {
                "name": "showyourwork",
                "tools": [
                    "showyourwork"
                ]
            }
        ],
        "latex-workshop.latex.tools": [
            {
                "name": "showyourwork",
                "command": "showyourwork",
                "args": [],
                "env": {}
            },
        ],
    }

Note that there should only be one set of outer braces.
In other words, remove the final outer brace in the first block above and the first outer brace in the second block above.


Debugging remote builds
-----------------------

Set ``CI=true`` environment variable to mimic the behavior on GitHub Actions locally.
More information on this soon.