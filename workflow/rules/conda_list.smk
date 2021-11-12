rule conda_list:
    """
    List the conda packages used in the current environment.

    """
    conda:
        posix(abspaths.user / "environment.yml")
    shell:
        "conda list"