import json
import jinja2
import subprocess
from pathlib import Path, PurePosixPath
from xml.etree.ElementTree import parse as ParseXMLTree
import re
import os
import numpy as np
import shutil


# TODO: Better way of getting abs paths?
WORKFLOW = Path(workflow.modules["showyourwork"].snakefile).absolute().parents[0]
USER = Path(workflow.basedir)


# Relative paths
TEX = Path("tex")
FIGURES = Path("figures")
TEMP = Path(".showyourwork")
GITHUB = Path(".github")


# Recognized figure extensions
FIGURE_EXTENSIONS = ["pdf", "png", "eps", "jpg", "jpeg", "gif", "svg", "tiff"]


# Dummy file dependency for figures w/ unknown parent scripts
UNKNOWN_SCRIPT = "unknown-script"


# Temporary tex files
TMPTEXFILE = ".showyourwork-xml-ms"
SYWTEXFILE = ".showyourwork-ms"


def posix(path):
    """
    Return the full POSIX path as a string.

    """
    return str(PurePosixPath(path))
