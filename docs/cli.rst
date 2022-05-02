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

.. |showyourwork_zenodo_help| raw:: html

    <span style="font-family:var(--pst-font-family-monospace);">showyourwork zenodo --help</span>


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


``showyourwork setup``
----------------------

.. admonition:: |showyourwork_setup_help|

    .. program-output:: showyourwork setup --help

The ``setup`` subcommand sets up an open source article repository from scratch
in the current working directory. This is an interactive command (unless you
provide the ``--yes`` option; see below). Let's step through what it does here.

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
or not the ``$ZENODO_TOKEN`` environment variable is set.


.. _syw_setup_step2a:

Step 2A
^^^^^^^

If the ``$ZENODO_TOKEN`` environment variable is not set,
you should see the following message:

.. raw:: html

    <pre>
    I didn't find a <span class="text-highlight">ZENODO_TOKEN</span> environment variable, so I'm not going to set up
    a Zenodo deposit for caching intermediate results. If you would like to enable
    this, please go to

        <a href="https://zenodo.org/account/settings/applications/tokens/new"><span class="text-highlight">https://zenodo.org/account/settings/applications/tokens/new</span></a>

    to create a new personal access token with deposit:actions and deposit:write
    scopes, store it in a local <span class="text-highlight">ZENODO_TOKEN</span> environment variable, and re-run this
    setup script.
    </pre>

    
The |showyourwork| workflow automatically caches the results of intermediate
steps in your pipeline on Zenodo, but only if it finds a ``$ZENODO_TOKEN`` 
environment
variable containing a valid Zenodo API token. If you would like to enable this
caching, exit out of the command by pressing ``Ctrl+C``, create a personal
access token on Zenodo with ``deposit:actions:`` and ``deposit:write`` scopes
at
`this page <https://zenodo.org/account/settings/applications/tokens/new>`_,
and save the token in an environment variable called ``$ZENODO_TOKEN``. Then 
re-run the ``setup`` command and check out :ref:`syw_setup_step2b` below.


.. _syw_setup_step2b:

Step 2B
^^^^^^^

If you set up a ``$ZENODO_TOKEN`` environment variable (see :ref:`syw_setup_step2a`), you should
instead see the following message:

.. raw:: html

    <pre>
    I found a <span class="text-highlight">ZENODO_TOKEN</span> environment variable, so I'm going to create a Zenodo
    deposit draft where intermediate results will be cached. In order for this to
    work on GitHub Actions, please go to

        <span class="text-highlight">https://github.com/rodluger/article/settings/secrets/actions/new</span>

    at this time and create a <span class="text-highlight">ZENODO_TOKEN</span> secret with your Zenodo access token.
    </pre>


As instructed in the message, go to your GitHub repository and create a "secret",
a secure variable that can be accessed by the GitHub Action that builds your
article on the cloud. Name this secret ``ZENODO_TOKEN`` and provide your
Zenodo API token (see above for details).

.. warning::

    Never commit your Zenodo API token (or any API token) directly to your
    repository!

You can read more about GitHub secrets (and the security measures
in place to prevent them from getting exposed to the outside world) at the 
`GitHub documentation <https://docs.github.com/en/actions/security-guides/encrypted-secrets>`_.

Press any key to bring up the next prompt. What you see next depends on whether
or not you specified the ``--overleaf`` option, and whether or not the environment
variables ``$OVERLEAF_EMAIL`` and ``$OVERLEAF_PASSWORD`` are set.


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

If you specified the ``--overleaf`` option (see :ref:`syw_setup_step3a`), but you
haven't configured your Overleaf credentials, you'll get the following message:

.. raw:: html

    <pre>
    It looks like you provided an Overleaf project id, but I didn't find an
    <span class="text-highlight">OVERLEAF_EMAIL</span> and/or an <span class="text-highlight">OVERLEAF_PASSWORD</span> environment variable, so I'm not
    going to set up Overleaf integration for this repository.
    </pre>

To allow |showyourwork| to push to/pull from your Overleaf project, create
the environment variables ``$OVERLEAF_EMAIL`` and ``$OVERLEAF_PASSWORD`` and 
populate them with your Overleaf email address and password, respectively;
then re-run the setup command.
Again, take care to never actually commit this information to your repository!


Step 3C
^^^^^^^

Finally, if you specified the ``--overleaf`` option *and* provided credentials
via the environment variables ``$OVERLEAF_EMAIL`` and ``$OVERLEAF_PASSWORD`` (see above), 
you'll get the following message:

.. raw:: html

    <pre>
    You provided an Overleaf project id, and I found both <span class="text-highlight">OVERLEAF_EMAIL</span> and
    <span class="text-highlight">OVERLEAF_PASSWORD</span> environment variables, so I'm going to set up Overleaf
    integration for this repository. In order for this to
    work on GitHub Actions, please go to

        <span class="text-highlight">https://github.com/rodluger/article/settings/secrets/actions/new</span>

    at this time and create <span class="text-highlight">OVERLEAF_EMAIL</span> and <span class="text-highlight">OVERLEAF_PASSWORD</span> secrets with
    your Overleaf credentials.
    </pre>

In order for the integration to work on GitHub Actions, you'll have to set the
repository secrets ``OVERLEAF_EMAIL`` and ``OVERLEAF_PASSWORD``, just as we
did for the ``ZENODO_TOKEN`` above.


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


Manual clean
^^^^^^^^^^^^

If ``showyourwork clean`` didn't remove all of the output, you can manually
delete all the programmatically-generated figures and datasets by removing
everything in the ``src/tex/figures`` and ``src/data`` folders
(assuming you're respecting the |showyourwork| conventions; see :doc:`layout`):

.. code-block:: bash

    rm -r src/tex/figures/**/*.*
    rm -r src/data/**/*.*

You may also have to manually remove the hidden ``.showyourwork`` folder, which
keeps track of repository metadata and caches certain files:

.. code-block:: bash

    rm -r .showyourwork


Deep clean
^^^^^^^^^^

If you want to start over from scratch, you can also delete the hidden ``.snakemake``
folder at the root of your repository:

.. code-block:: bash

    rm -r .snakemake

This houses the ``conda`` environments for your build (among other things), so
deleting it will force a re-install of all packages used in your workflow.

Finally, there's one more hidden folder to know about, a ``.showyourwork``
folder located in your ``$HOME`` path, which also houses ``conda`` environments
used at different stages of the build step. You can safely remove it at any time
(at the cost of a longer runtime the next time you execute ``showyourwork``):

.. code-block:: bash

    rm -r ~/.showyourwork


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


.. _syw_zenodo:

``showyourwork zenodo``
-----------------------

.. admonition:: |showyourwork_zenodo_help|

    .. program-output:: showyourwork zenodo --help

Utilities for creating, deleting, and publishing the Zenodo deposit drafts used
to cache intermediate results from your workflow; see :doc:`zenodo`.