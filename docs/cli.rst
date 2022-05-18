.. |showyourwork_help| raw:: html

    <span style="font-family:var(--pst-font-family-monospace);">showyourwork --help</span>

.. |showyourwork_setup_help| raw:: html

    <span style="font-family:var(--pst-font-family-monospace);">showyourwork setup --help</span>

.. |showyourwork_build_help| raw:: html

    <span style="font-family:var(--pst-font-family-monospace);">showyourwork build --help</span>

.. |showyourwork_clean_help| raw:: html

    <span style="font-family:var(--pst-font-family-monospace);">showyourwork clean --help</span>

.. |showyourwork_tarball_help| raw:: html

    <span style="font-family:var(--pst-font-family-monospace);">showyourwork tarball --help</span>

.. |showyourwork_cache_help| raw:: html

    <span style="font-family:var(--pst-font-family-monospace);">showyourwork cache --help</span>
 

Command line interface
======================

The |showyourwork| package implements a single command-line utility:
``showyourwork``, which allows users to set up, configure, and build their
open source article. Below we describe this command and discuss its various
subcommands.


``showyourwork``
----------------

.. admonition:: |showyourwork_help|

    .. program-output:: showyourwork --help

Running |showyourwork| (without any arguments) is a shortcut for ``showyourwork build``
(see :ref:`syw_build` below).


.. _syw_setup:

``showyourwork setup``
----------------------

.. admonition:: |showyourwork_setup_help|

    .. program-output:: showyourwork setup --help

The ``setup`` subcommand sets up an open source article repository from scratch
in the current working directory. This is an interactive command (unless you
provide the ``--yes`` or ``--quiet`` options; see below). 
Let's step through what it does here.

To set up a new open source article repository, run

.. raw:: html

    <pre>
    showyourwork setup <span class="text-highlight">$USER/$REPO</span>
    </pre>

where you should replace ``$USER`` with your GitHub user name and 
``$REPO`` with the name of your new article repository. For definiteness,
here we'll use my user name (``rodluger``) and we'll call our repository
``article``.


Step 1
^^^^^^

Running the ``setup`` command as above should bring up the following prompt:

.. raw:: html

    <pre>
    Let's get you set up with a new repository. I'm going to create a folder called

        <span class="text-highlight">article</span>

    in the current working directory. If you haven't done this yet, please visit

        <a href="https://github.com/new"><span class="text-highlight">https://github.com/new</span></a>

    at this time and create an empty repository called

        <span class="text-highlight">rodluger/article</span>
    </pre>


As requested, if you haven't yet created the remote repository, go to
`github.com/new <https://github.com/new>`_ in your browser to create an empty 
repository of the same name. There's no need to create a README, gitignore file, 
or LICENSE at this time, as |showyourwork| will set those up for you.

Press any key to bring up the next prompt. What you see next depends on whether
or not you provided the ``--cache`` option.


.. _syw_setup_step2a:

Step 2A
^^^^^^^

If you didn't request Zenodo caching functionality (see below), you'll see the
following the message:

.. raw:: html

    <pre>
    You didn't provide a caching service (via the <span class="text-highlight">--cache</span>
    command-line option), so I'm not going to set up remote caching for this repository.
    </pre>


.. _syw_setup_step2b:

Step 2B
^^^^^^^

If instead you passed the ``--cache`` flag, you'll see the following message:

.. raw:: html

    <pre>
    You requested remote caching on Zenodo, so I'm going to create a deposit draft where 
    intermediate results will be cached. Please make sure at this time that you have 
    defined the <span class="text-highlight">ZENODO_TOKEN</span> environment variable containing your API key for Zenodo. 
    If you don't have one, visit

        <span class="text-highlight">https://zenodo.org/account/settings/applications/tokens/new</span>

    to create a new personal access token with <span class="text-highlight">deposit:actions</span> and <span class="text-highlight">deposit:write</span> 
    scopes and store it in the environment variable <span class="text-highlight">ZENODO_TOKEN</span>. In order for 
    this to work on GitHub Actions, you'll also have to visit

        <span class="text-highlight">https://github.com/tmp/tmp/settings/secrets/actions/new</span>

    at this time to create a <span class="text-highlight">ZENODO_TOKEN</span> secret with your API access token.
    </pre>

Note that, in addition to the ``--cache`` flag, which enables caching on Zenodo,
users may also provide the ``--sandbox`` flag, which switches the host to Zenodo
Sandbox. Zenodo Sandbox behaves in exactly the same way as Zenodo, but it is
explicitly meant as a test bed for dataset archiving. While deposits on Sandbox get assigned
DOIs, they are no *actual* registered DOIs and have a limited lifespan.
Sandbox is therefore a great choice for debugging and development; read more about
it at :doc:`zenodo`. Note that if you choose the ``--sandbox`` option, you'll need
a Zenodo Sandbox API token stored int the ``SANDBOX_TOKEN`` environment variable
and GitHub Actions secret.

.. warning::

    Never commit your Zenodo API token (or any API token) directly to your
    repository!

You can read more about GitHub secrets (and the security measures
in place to prevent them from getting exposed to the outside world) at the 
`GitHub documentation <https://docs.github.com/en/actions/security-guides/encrypted-secrets>`_.

Press any key to bring up the next prompt. What you see next depends on whether
or not you specified the ``--overleaf`` option.


.. _syw_setup_step3a:

Step 3A
^^^^^^^

If you didn't pass the ``--overleaf`` option, you'll see the following:

.. raw:: html

    <pre>
    You didn't provide an Overleaf project id (via the <span class="text-highlight">--overleaf</span> command-line
    option), so I'm not going to set up Overleaf integration for this repository.
    </pre>

If you would like to set up integration with an Overleaf project (see :doc:`overleaf`),
hit ``Ctrl+C`` and run

.. code-block:: bash

    showyourwork setup --overleaf=62150dd16134ef045f81d1c8

where you should replace ``62150dd16134ef045f81d1c8`` with the 24-character id 
of a new (blank) Overleaf project. Once you create a new Overleaf project, you
can grab the id from the last bit of the project's URL. Note that |showyourwork|
requires the Overleaf project to be empty, otherwise it will refuse to set up
the integration. For more information on how this integration works, and what
to do if you have an existing Overleaf project you'd like to integrate with
|showyourwork|, please see :doc:`overleaf`.


Step 3B
^^^^^^^

If you specified the ``--overleaf`` option (see :ref:`syw_setup_step3a`),
you'll get the following message:

.. raw:: html

    <pre>
    You provided an Overleaf project id, so I'm going to set up Overleaf integration 
    for this repository. Please make sure at this time that you have defined the 
    <span class="text-highlight">OVERLEAF_EMAIL</span> and <span class="text-highlight">OVERLEAF_PASSWORD</span> environment variables. In order for this to
    work on GitHub Actions, please go to

        <span class="text-highlight">https://github.com/tmp/tmp/settings/secrets/actions/new</span>

    at this time and create <span class="text-highlight">OVERLEAF_EMAIL</span> and <span class="text-highlight">OVERLEAF_PASSWORD</span> secrets with your 
    Overleaf credentials.
    </pre>

To allow |showyourwork| to push to/pull from your Overleaf project, create
the environment variables ``$OVERLEAF_EMAIL`` and ``$OVERLEAF_PASSWORD`` and 
populate them with your Overleaf email address and password, respectively;
then re-run the setup command.
Again, take care to never actually commit this information to your repository!


Step 4
^^^^^^

Finally, press any key to generate the repository. This will create a new folder
in the current working directory with the same name as your repo (``article``, in
the example above) and set up ``git`` tracking for it. Note that the first time
you commit and push your changes to the GitHub repository, you'll have to specify
the upstream branch as follows:

.. code-block:: bash

    git push -u origin main


.. _syw_build:

``showyourwork build``
----------------------

.. admonition:: |showyourwork_build_help|

    .. program-output:: showyourwork build --help

Run this command to build the article in the current working directory. Note that
you must run this command from the top level of the repository (an error will
be thrown otherwise). The command accepts any number of arguments, all of which
are forwarded to ``snakemake``. 
By default, ``showyourwork`` passes the following arguments to ``snakemake``:

.. code-block:: bash

    -c1 --use-conda --reason --cache

Some of these, like the number of cores, can be overridden. For example, you
may run

.. code-block:: bash

    showyourwork build -c2

to run the workflow using two cores (see the `snakemake docs <https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_
for details). Additional arguments can also be provided, like ``--verbose`` to increase
the verbosity of the Snakemake logs (see :doc:`logging`), or ``--force`` and ``--forceall`` to
force the re-execution of the rule that builds the manuscript or *all* of the rules
in the workflow, respectively (regardless of whether the outputs are up to date
or not). Positional arguments are also allowed; for instance, to only build a specific
figure, you may run, e.g.,

.. code-block:: bash

    showyourwork build --force src/tex/figures/figure.pdf

You can check out the complete list of Snakemake arguments and options
at the `snakemake documentation <https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options>`_.

.. warning::

    Not all Snakemake options are compatible with |showyourwork|. If you
    run into issues when specifying custom options, please 
    `let us know <https://github.com/showyourwork/showyourwork/issues/new>`_.

Note that the build process in |showyourwork| happens in two steps, each of
which executes a separate Snakemake workflow. The first
step is a preprocessing step that parses the user config file and does a quick
first-pass compiling of the TeX manuscript to look for ``\includegraphics``
and ``\script`` calls, which it uses to build the graph of dependencies for
your article. The second step is the main step, in which all of the dependencies
are built (if needed) and the final article PDF is generated. Arguments
passed to ``showyourwork build`` are ingested *only* during the second step.

Finally, |showyourwork| takes full advantage
of the dependency tracking and caching functionality of Snakemake. When
running ``showyourwork build``, only files whose upstream dependencies have
changed (since the last build) will be re-generated. This is true *even when
running on GitHub Actions*; the ``showyourwork-action`` caches results across
runs to minimize compute time for the build. We even go a step further, and
extend the Snakemake functionality to allow caching of intermediate
dependencies on Zenodo; read about it at :doc:`zenodo`.


.. _syw_clean:

``showyourwork clean``
----------------------

.. admonition:: |showyourwork_clean_help|

    .. program-output:: showyourwork clean --help

This command removes all of the output from previous ``showyourwork build``
steps. Depending on the state of your repository, and if there are errors in
your config file or missing dependencies, this command may fail silently, in
which case some of the output may remain after running it.


Force clean
^^^^^^^^^^^

If ``showyourwork clean`` didn't remove all of the output, you can force
the deletion of all the programmatically-generated figures and datasets by 
passing the ``--force`` option, which will remove
everything in the ``src/tex/figures``, ``src/data``, and temporary
``.showyourwork`` folders
(assuming you're respecting the |showyourwork| conventions; see :doc:`layout`).


Deep clean
^^^^^^^^^^

If you want to start over from scratch, you can also pass the ``--deep`` option
for a deep clean. This removes
the hidden ``.snakemake`` folder, which houses the ``conda`` environments 
for your build (among other things), so
deleting it will force a re-install of all packages used in your workflow.
This option also removes the ``.showyourwork``
folder located in your ``$HOME`` path, which also houses ``conda`` environments
used at different stages of the build step. You can safely remove it at any time
(at the cost of a longer runtime the next time you execute ``showyourwork``).


.. _syw_tarball:

``showyourwork tarball``
------------------------

.. admonition:: |showyourwork_tarball_help|

    .. program-output:: showyourwork tarball --help

Like ``build``, the ``showyourwork tarball`` command builds your article, but
also gathers all of the relevant files needed to build it using a standard
TeX engine into a tarball called ``arxiv.tar.gz``. It's named that because
you should be able to directly upload this tarball when submitting a paper
to the `arXiv <https://arxiv.org/>`_ article service.


.. _syw_cache:

``showyourwork cache``
----------------------

.. admonition:: |showyourwork_cache_help|

    .. program-output:: showyourwork cache --help

    .. raw:: html

        <br/>

    **Subcommand documentation:**

    .. program-output:: showyourwork cache create --help

    .. program-output:: showyourwork cache delete --help

    .. program-output:: showyourwork cache freeze --help

    .. program-output:: showyourwork cache publish --help

Utilities for creating, deleting, and publishing the Zenodo deposit drafts used
to cache intermediate results from your workflow; see :doc:`zenodo`.