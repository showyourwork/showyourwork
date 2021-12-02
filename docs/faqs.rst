FAQs
====

Below is a non-exhaustive list of frequently asked questions. If you don't
see what you're looking for here, check out the `issues page <https://github.com/rodluger/showyourwork/issues>`_ or feel
free to `suggest your own FAQ <https://github.com/rodluger/showyourwork/edit/main/docs/faqs.rst>`_.

What if I already have a repository?
------------------------------------

We recommend starting your ``showyourwork`` project from scratch by
instantiating the `repository template <https://github.com/rodluger/showyourwork-template/generate>`_.
If you already have a project in progress (or even completed), you'll have to do
a bit of re-structuring to get it to conform to the ``showyourwork`` layout. This means
placing everything inside a top-level ``src`` folder, and within that, figure scripts
within the ``figures`` folder, static files within the ``static`` folder, and the TeX file
directly under ``src``. Next, you'll have to add the ``showyourwork`` submodule by running

.. code-block:: bash

    git submodule add https://github.com/rodluger/showyourwork

Finally, copy the files ``Snakefile``, ``Makefile``, ``environment.yml``, ``showyourwork.yml``, 
and ``.github/workflow/showyourwork.yml``
from the `showyourwork template <https://github.com/rodluger/showyourwork-template>`_,
and edit them as needed. See :doc:`the layout guide <layout>` for details.


What if I'm using Overleaf?
---------------------------

We currently don't support Overleaf projects. We're thinking about how to do this properly, so
please check back soon for more information. In the meantime, you can develop both your ``showyourwork``
repository and your Overleaf project simultaneously, copying the tex file over to GitHub and the figure
output over to Overleaf manually as needed.


What if I don't use Python?
---------------------------

By default, ``showyourwork`` expects figure scripts to be Python ``.py`` scripts.
Other languages are also supported, but users must provide explicit instructions
on how to produce figures from the corresponding figure scripts. See the :doc:`custom`
page for details on how to get ``showyourwork`` to recognize non-Python scripts.


What if I don't use ``includegraphics`` calls?
----------------------------------------------

``showyourwork`` inspects calls to ``includegraphics`` and ``label`` within ``figure``
environments to infer figure dependencies in the workflow. However,
many people use commands like ``plotone`` and ``plottwo`` or user-defined shorthand for
including figures in LaTeX. **This should still work!** As long as these shortcuts
call ``includegraphics`` somewhere along the line (which they do), ``showyourwork``
should work just fine. Just remember to always label your figures!


Can I nest figure scripts inside folders?
-----------------------------------------

The suggested workflow is to place all your figures directly under the ``src/figures``
directory. However, you may nest them under additional subfolders in that directory.
Note that you must change your figure labels accordingly (i.e., if your figure script is 
``src/figures/some_folder/script.py``, you should label the corresponding figure in LaTeX
as ``fig:some_folder/script``.) The other catch is that your script will be run from
the ``src/figures`` directory, so just make sure the file is saved to the current
working directory (not the directory the script is in). This is the default behavior
when doing I/O in Python, so it shouldn't generally be a problem.


What if I'm on Windows?
-----------------------

I don't know! Does it work? Please give me your feedback `here <https://github.com/rodluger/showyourwork/issues/33>`_.
I'm hoping to add unit tests for Windows and add explicit support for it soon,
but any help is appreciated!


How do I debug a run?
---------------------

Sometimes it can be difficult to figure out why a run failed from the terminal
output. Try increasing the verbosity of the output to obtain more informative 
error messages by specifying ``verbose: true`` in your ``showyourwork.yml``
file. If you've fixed an issue and are still getting the same error, you might try
deleting the ``.showyourwork`` folder at the top level of the repo to get rid
of the build cache. You can also try running ``make clean``, although that will
delete all your build output.

If your local build is passing but the remote (GitHub Actions) build
is failing, this could be due to missing files on the remote. This can 
happen if you forgot to ``git add`` and push one or more 
files to the remote repository. Some files are not tracked by
default (like things in the ``data`` directory, non-Python scripts in the
``figures`` directory, or anything other than ``tex`` and ``bib`` files
in the ``src`` directory), so you'll have to ``git add -f`` (force add) them.

When a remote build fails, an artifact is uploaded called ``showyourwork-output``
(see `the docs <https://docs.github.com/en/actions/managing-workflow-runs/downloading-workflow-artifacts>`_ 
for details on how to download it). This is a tarball containing all files in the
root of your repository and all files directly under the ``src`` directory
at the time the build failed, excluding large (> 5MB) files. Poking around here
can sometimes help identify the cause of the failure. You can also try
increasing the verbosity in the ``.github/workflows/showyourwork.yml`` file
(see :ref:`github_action_verbose`.)

Finally, one thing you should routinely do is :doc:`update showyourwork <faqs>`.
We're constantly fixing issues on your side, so doing this might just fix your problem!

If none of this helps, please check out the 
`issues <https://github.com/rodluger/showyourwork/issues?q=is%3Aissue>`_
page for ``showyourwork`` to see if someone has run into the same problem before.
Feel free to open a new one if you don't see what you're looking for.


What if I don't use AASTeX?
---------------------------

If you're not using the `AASTeX <https://journals.aas.org/aastexguide/>`_ or 
some of the other built-in templates (see :doc:`Custom workflows <custom>`), 
things should still work as long as you include the
relevant class and style files in the ``src/`` directory alongside your texfile.
You won't get the snazzy ``showyourwork`` logo at the top of the page, but
everything else should still work. We're planning on adding explicit support for
other templates, so please check back soon for more or open 
`an issue <https://github.com/rodluger/showyourwork/issues?q=is%3Aissue>`_.


I get a warning saying the Zenodo upload failed.
------------------------------------------------

If you don't have the right authentication, and a workflow attempts to 
publish a deposit to Zenodo under a certain ``id``, you will get a warning
saying something along the lines of 

.. code-block:: bash

    Error: Unable to upload <file-name> to Zenodo.

and

.. code-block:: bash

    Zenodo error 401: The server could not verify that you are authorized to access the URL requested. You either supplied the
    wrong credentials (e.g. a bad password), or your browser doesn't understand how to supply the credentials required.

This can happen if you forgot to set your Zenodo API token environment variable
(see the :ref:`token_name <zenodo.dataset.token_name>` config setting for details)
or if you've cloned a third-party repository and are trying to reproduce their 
results locally. In the latter case, the easiest workaround is to run

.. code-block:: bash

    make fast

which will skip the generation & upload step for any file that can instead be
downloaded from Zenodo. Alternatively, you can change the relevant ``id``s in the
``shoyourwork.yml`` config file to *version* ids, which correspond to static
(download-only) entries (:ref:`read more about that here <zenodo.dataset.id>`),
or change them to concept ids that you have access to (you can obtain one
by running ``make reserve``).