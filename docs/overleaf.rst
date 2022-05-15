Overleaf integration
====================

.. warning::

    Overleaf integration is still an **experimental** feature.
    The user interface may change in future releases.
    If you run into problems or have any suggestions to improve
    this functionality, please 
    `open an issue on GitHub <https://github.com/showyourwork/showyourwork/issues/new>`__.

Overview
--------

|showyourwork| now allows users to integrate their projects 
with `Overleaf <https://www.overleaf.com>`__, which can greatly facilitate collaborative
article writing in LaTeX.
While Overleaf supports integration with git and GitHub in 
`a few different ways <https://www.overleaf.com/learn/how-to/Using_Git_and_GitHub>`__,
none of these are *quite* right for an Overleaf-|showyourwork| bridge (you
can read all about why `here <https://github.com/showyourwork/showyourwork/issues/22>`__).

Instead, whenever ``showyourwork`` is executed (either locally or on GitHub Actions),
a specific set of files is synced between the git repository and the Overleaf project.
The user can choose which files get synced, but the default is to push figure output
to Overleaf and pull TeX files from Overleaf.
In order to prevent merge conflicts and data loss, syncing only happens in one direction:
by default,
figures only ever get **pushed to** Overleaf, while TeX files only ever get **pulled from**
Overleaf to the git repository.
This means that if users wish to integrate their project with Overleaf, they must only
make changes to the TeX file(s) on Overleaf. Similarly, users should not manually change/edit/upload
figures to Overleaf, and instead let |showyourwork| do the grunt work.

In the event of a merge conflict---for example, if the user edited ``ms.tex`` locally---the 
|showyourwork| build will acknowledge that and fail, refusing to overwrite
the local changes. See :ref:`conflicts` below for details.

Below, we'll discuss how to set up the integration with Overleaf and take a look at how
to configure the syncing.


Setup
-----

Let's set up Overleaf integration for a new repository. It is much, much, much
easier to set up Overleaf integration for new projects, so here we'll create both a 
new GitHub repository **and** a new Overleaf project for our article.
If you want to set up Overleaf integration for an existing project, see :ref:`existing`.

To start, create a new repository on GitHub by visiting `github.com/new <https://github.com/new>`__. 
For definiteness,
in the example below we'll create the repository ``article`` under my user
name (``rodluger``).

Next, create a new project on `Overleaf <https://www.overleaf.com/project>`__ from
the **Blank Project** template. You'll be directed to the new project, whose URL
will look something like

.. code-block:: text

    https://www.overleaf.com/project/6272c02ffe09ce2c9a5f0ff6

That last bit of the URL (``6272c02ffe09ce2c9a5f0ff6``) is the 24-character
project ID, which uniquely identifies your project; we'll use that in just
a moment. On your project page, you should see a single file called ``main.tex``
containing some minimal TeX stuff. Don't edit anything in it, as it will get
overwritten momentarily by |showyourwork|.

Next, on your machine, define the environment variables ``$OVERLEAF_EMAIL``
and ``$OVERLEAF_PASSWORD`` with your Overleaf credentials. Unfortunately,
Overleaf does not yet support token-based authentication, so the only way
to grant |showyourwork| access to your project is by providing the email
and password for your Overleaf account. As with the other credentials required
by |showyourwork|, you'll also need to create corresponding
GitHub Actions secrets at the following URL 

.. code-block:: text

    https://github.com/rodluger/article/settings/secrets/actions/new

where you should replace ``rodluger/article`` with your GitHub username
and the name of the repository.

.. warning::

    Never commit your Overleaf credentials (or any credentials) directly to your
    repository!

Finally, on your local machine, set up a new |showyourwork| repository by running

.. code-block:: bash

    showyourwork setup rodluger/article --overleaf=6272c02ffe09ce2c9a5f0ff6

where, again, you should replace ``rodluger/article`` with your repository slug
and ``6272c02ffe09ce2c9a5f0ff6`` with your Overleaf project ID. You should see
a couple messages saying the Overleaf repository is being set up. Once that's
done, if you navigate to your Overleaf project in the browser, you'll see that
the single TeX file ``main.tex`` is gone and has been replaced by the following
files:

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
          <li class="folder">figures
            <ul>
                <li>.gitignore</li>
            </ul>
          </li>
          <li>.gitignore</li>
          <li>bib.bib</li>
          <li>ms.tex</li>
          <li>showyourwork.sty</li>
        </ul>
      </div>

These, in fact, are the same files as in the ``src/tex`` folder of your repository
(see :doc:`layout`); |showyourwork| will keep your Overleaf project up to date
with the contents of that folder (more on this below). 
Note that the TeX manuscript is now called 
``ms.tex`` (the default name in |showyourwork|).

Returning to our local |showyourwork| repository, if you open the config file
``showyourwork.yml``, you'll see that the ``setup`` command populated the
``overleaf`` field with some stuff:

.. code-block:: yaml

    overleaf:
        id: 6272c02ffe09ce2c9a5f0ff6
        push: 
            - src/tex/figures
        pull:
            - src/tex/ms.tex
            - src/tex/bib.bib

In addition to your Overleaf project ID, it has also defined some files and folders
in the ``push`` and ``pull`` fields. To understand what these mean, read on!


Pushing and pulling
-------------------

Coming soon! Only enabled for the ``main`` branch.


.. _conflicts:

Managing conflicts
------------------

Coming soon!



.. _existing:

Integrating existing projects
-----------------------------

Coming soon!