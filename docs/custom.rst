Custom workflows
================

By default, the workflow defined in the ``Snakefile`` looks like this:

.. code-block:: python

    # Import the showyourwork module
    module showyourwork:
    snakefile:
        "showyourwork/workflow/Snakefile"
    config:
        config


    # Use all default rules
    use rule * from showyourwork


Custom dependencies
-------------------

Download a dataset and make it a dependency of a particular figure:

.. code-block:: python

    # Custom rule to download a dataset
    rule my_dataset:
        output:
            report("src/figures/my_dataset.dat", category="Dataset")
        shell:
            "curl https://zenodo.org/record/5187276/files/fibonacci.dat --output {output[0]}"

    # Subclass the `figure` rule to specify that `src/data/dataset.dat`
    # is a dependency of `src/figures/my_figure.pdf`
    use rule figure from showyourwork as my_figure with:
        input:
            "src/figures/my_figure.py",
            "src/figures/my_dataset.dat",
            "environment.yml"
        output:
            report("src/figures/my_figure.pdf", category="Figure")


Custom figure scripts
---------------------

Specify a custom script for a figure. Useful when ``showyourwork`` can't
automatically determine the figure script, such as when a figure is
included outside of a ``figure`` environment:

.. code-block:: python

    # Subclass the `figure` rule to specify that the figure
    # `src/figures/custom_figure.pdf` is generated from the script
    # `src/figures/custom_script.py`
    use rule figure from showyourwork as custom_figure with:
        input:
            "src/figures/custom_script.py",
            "environment.yml"
        output:
            report("src/figures/custom_figure.pdf", category="Figure")


Override the internal ``figure`` rule completely:

.. code-block:: python

    rule custom_figure:
        input:
            "src/figures/custom_script.py",
            "environment.yml",
        output:
            report("src/figures/custom_figure.pdf", category="Figure")
        conda:
            "environment.yml"
        bash:
            "cd src/figures && python custom_script.py"
