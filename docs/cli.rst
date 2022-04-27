.. raw:: html

    <style>
        .text-highlight{
            font-weight:bold; 
            color:#204a87;
        }
    </style>

.. |showyourwork_help| raw:: html

    <span style="font-family:var(--pst-font-family-monospace);">showyourwork --help</span>

.. |showyourwork_setup_help| raw:: html

    <span style="font-family:var(--pst-font-family-monospace);">showyourwork setup --help</span>


Command line interface
======================

The ``showyourwork`` package implements a single command-line command:
``showyourwork``, which allows users to set up, configure, and build their
open source article. Below we list the various subcommands and what they do.


``showyourwork``
----------------

.. admonition:: |showyourwork_help|

    .. program-output:: showyourwork --help

Running ``showyourwork`` (without any arguments) is a shortcut for ``showyourwork build``
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
or LICENSE at this time, as ``showyourwork`` will set those up for you.

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

    
The ``showyourwork`` workflow automatically caches the results of intermediate
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
can grab the id from the last bit of the project's URL. Note that ``showyourwork``
requires the Overleaf project to be empty, otherwise it will refuse to set up
the integration. For more information on how this integration works, and what
to do if you have an existing Overleaf project you'd like to integrate with
``showyourwork``, please see :doc:`overleaf`.


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

To allow ``showyourwork`` to push to/pull from your Overleaf project, create
the environment variables ``$OVERLEAF_EMAIL`` and ``$OVERLEAF_PASSWORD`` and 
populate them with your Overleaf email address and password, respectively;
then re-run the setup command.
Again, take care to never actually commit this information to your repository!


Step 3C
^^^^^^^

If you followed the instructions above, you'll see the following message:

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
the example above) and set up ``git`` tracking for it.


.. _syw_build:

``showyourwork build``
----------------------

.. note::

    Coming soon!