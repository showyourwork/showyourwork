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

Currently, ``showyourwork`` requires figure scripts to be Python ``.py`` scripts.
Support for other languages -- such as ``julia``, ``Jupyter`` notebooks, or even
just plain ``bash`` is coming soon! In the meantime, check out :doc:`custom` for
information on how to override this behavior manually to specify your own figure
generation rules.


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


How do I debug a run?
---------------------

Sometimes it can be difficult to figure out why a run failed from the terminal
output, especially if it's an error during the LaTeX build. Try increasing the
verbosity of the output to obtain more informative error messages:

.. code-block:: bash

    snakemake -c1 --use-conda ms.pdf --config verbose=true

If that doesn't help, check out the `issues <https://github.com/rodluger/showyourwork/issues?q=is%3Aissue>`_
page for ``showyourwork`` to see if someone has run into the same problem before.
Feel free to open a new one if you don't see what you're looking for.


What if I don't use AASTeX?
---------------------------

If you're not using the `AASTeX <https://journals.aas.org/aastexguide/>`_ LaTeX 
template for AAS journals, things should still work as long as you include the
relevant class and style files in the ``src/`` directory alongside your texfile.
You won't get the snazzy ``showyourwork`` logo at the top of the page, but
everything else should still work. We're planning on adding explicit support for
other templates (like MNRAS and A&A), so please check back soon for more.


Why do I keep getting a Bad Gateway error?
------------------------------------------

On GitHub Actions you may occasionally run into a ``502: Bad Gateway`` HTTP
error. These are sporadic and temporary, and are usually due to issues
with downloading the (very large) TeX distributions from ``archive.org``. 
To commiserate with others experiencing this issue, `see here <https://github.com/tectonic-typesetting/tectonic/issues/765>`_.
Other than that, the best thing to do is simply wait a bit and re-run the failed job!