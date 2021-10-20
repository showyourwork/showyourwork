"""
Sets global constants for the workflow.

"""
from pathlib import Path, PurePosixPath
import os
from sphinx_mock import *


_workflow = workflow


def posix(path):
    """Convert a path object to a string containing the posix path."""
    return str(PurePosixPath(path))


class abspaths:
    """Absolute paths used throughout the workflow."""

    #: Path to the showyourwork workflow
    workflow = Path(_workflow.modules["showyourwork"].snakefile).absolute().parents[0]

    #: Path to the user's repository root
    user = Path(_workflow.basedir).absolute()


class relpaths:
    """Relative paths used throughout the workflow."""

    #: Path to the showyourwork workflow
    tex = Path("src")

    #: Path to the showyourwork workflow
    figures = Path("src") / "figures"

    #: Path to the showyourwork workflow
    dot_github = Path(".github")

    #: Path to the showyourwork workflow
    temp = Path(".showyourwork")
    if not temp.exists():
        os.mkdir(str(temp))