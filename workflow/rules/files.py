"""
Miscellaneous files used throughout the workflow.

"""


class files:
    """Houses lists of support files used throughout the workflow."""

    #: The `.zenodo` files in which we store the URL for the datasets
    #: We make the PDF generation depend on these so we always have
    #: access to the URLs for the icon generation.
    #: These are updated in `zenodo.smk`.
    dot_zenodo = []

    #: Auxiliary files we copy over to the user's `src/` directory.
    aux = [
        posix(relpaths.tex / file.name)
        for file in (abspaths.workflow / "resources" / "tex").glob("*.*")
    ]

    #: Class-specific auxiliary files.
    cls = {}
    for folder in (abspaths.workflow / "resources" / "classes").glob("*"):
        cls[folder.name] = [
            file.name
            for file in (abspaths.workflow / "resources" / "classes" / folder).glob(
                "*.*"
            )
        ]

    #: Dummy file dependency for figures w/ unknown parent scripts.
    unknown = "unknown-script"

    #: Temporary tex file.
    tmp_xml = ".showyourwork-xml-ms"

    #: Temporary tex file.
    tmp_syw = ".showyourwork-ms"

    #: Tectonic command
    if config["tectonic_latest"]:
        tectonic = [posix(relpaths.temp / "tectonic")]
    else:
        tectonic = []

    #: Figures that are allowed directly in the ``src/`` directory
    special_figures = ["orcid-id.png", "showyourwork.pdf"]
