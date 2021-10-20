"""
Miscellaneous files used throughout the workflow.

"""
from sphinx_mock import *


class files:
    """Houses lists of support files used throughout the workflow."""

    #: The manually-uploaded datasets for this article that must be downloaded
    zenodo_files_manual = []

    #: The showyourwork-managed datasets for this article
    zenodo_files_auto = []

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

    if config["tectonic_latest"]:
        #: Tectonic command
        tectonic = [posix(relpaths.temp / "tectonic")]
    else:
        #: Tectonic command
        tectonic = []

    #: Figures that are allowed directly in the ``src/`` directory
    special_figures = ["orcid-id.png", "showyourwork.pdf"]
