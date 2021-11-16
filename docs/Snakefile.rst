The ``Snakefile``
=================

This is the file defining the ``Snakemake`` workflow for the article
build. You can read more about these kinds of files in the
`Snakemake documentation <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html>`_.
The syntax used in this file is a superset of ``python``, meaning it supports
almost all the usual ``python`` syntax, plus several new directives.

Let's take a look at the default contents of the ``Snakefile``:

.. code-block:: python

    # User config
    configfile: "showyourwork.yml"


    # Import the showyourwork module
    module showyourwork:
        snakefile:
            "showyourwork/workflow/Snakefile"
        config:
            config


    # Use all default rules
    use rule * from showyourwork

There are several things happening here. First, this line

.. code-block:: python

    # User config
    configfile: "showyourwork.yml"

tells ``Snakemake`` that the configuration file for the workflow is
``showyourwork.yml`` (see `the showyourwork.yml file <config.html>`_).
Next, 

.. code-block:: python

    # Import the showyourwork module
    module showyourwork:
        snakefile:
            "showyourwork/workflow/Snakefile"
        config:
            config

imports the ``showyourwork`` workflow, and 

.. code-block:: python

    # Use all default rules
    use rule * from showyourwork

tells ``Snakemake`` to use all of the rules defined in that workflow (this
is similar to a ``from showyourwork import *`` command in vanilla ``python``).

You can edit this file to customize your workflow, such as to
define custom ``rule`` directives to generate figures (see the 
`custom figure scripts example <custom.html#custom-figure-scripts>`_)
or to read a programmatically-generated config file
(see the 
`many dependencies example <custom.html#many-many-dependencies>`_).

You can do more advanced stuff, too, but that likely requires an
understanding of the ``showyourwork`` workflow and what all of its
rules do. If you're interested, have a look at the
`Developer API <api.html>`_ page; in principle, all of the rules
defined there can be overridden in your workflow.