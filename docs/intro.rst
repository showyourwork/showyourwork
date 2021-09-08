Introduction
============

The ``showyourwork`` workflow is intended to help authors publish the code that generated
the figures and results in a scientific article. It ensures that the compiled article
PDF is always in sync with all of the code used to generate it. It does this
automatically—and seamlessly—with the help of the `Snakemake workflow management system <https://snakemake.readthedocs.io>`_,
the `tectonic typesetting engine <https://tectonic-typesetting.github.io>`_, and
`Github Actions CI <https://github.com/features/actions>`_.

.. raw:: html

    <div class="admonition note">
        <p class="admonition-title">The showyourwork philosophy</p>
        <p>Scientific papers should exist as GitHub repositories comprised of
        LaTeX files, figure scripts, rules to access datasets, a platform/environment specification,
        <span style="font-weight:bold;">and nothing else</span>.
        Anyone should be able to re-generate the article PDF from scratch at
        the click of a button.
        </p>
    </div>

Within the ``showyourwork`` workflow, scientific articles exist as GitHub repositories
with a :doc:`specific layout <layout>`. Whenever new commits are pushed to the remote
repository, a GitHub action is triggered that automatically builds the article from the
input figure scripts, manuscript files, and
`conda environment file <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_,
following the instructions specified in the `Snakefile <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html>`_:

.. raw:: html

    <div align="center" style="margin-bottom: 17.25px;">
        <img src="https://raw.githubusercontent.com/rodluger/showyourwork/img/workflow.png" width="80%"/>
    </div>

The article PDF (along with a tarball containing all of the output) is then pushed to a special branch
(usually ``main-pdf``) on the repository. This article is decorated with badges linking to the exact
versions of the files on GitHub used to generate it.

Thanks to the magic of ``Snakemake``, ``showyourwork`` is both lightweight—it should work out-of-the-box for most users—and highly
customizable. It also uses intelligent caching to never re-run things it doesn't have to (like figure scripts that haven't changed),
even when running on the cloud.

To get started with ``showyourwork``, check out the :doc:`quickstart tutorial <quickstart>`.
You should also read about
the :doc:`showyourwork GitHub action <action>`, how to :doc:`build your article locally <local>`,
and how to :doc:`customize your workflow <custom>`.

You should also spend some time browsing through the :doc:`FAQs page <faqs>`. Since ``showyourwork``
is itself a work in progress, new features are still being added frequently. If you spot a bug,
have a question, or would like ``showyourwork`` to do something it doesn't currently support,
please feel free to raise a `GitHub issue <https://github.com/rodluger/showyourwork/issues>`_.
