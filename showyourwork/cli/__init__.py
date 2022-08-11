"""
Command line interface for ``showyourwork``.

This module provides the command line interface for ``showyourwork``,
consisting of the ``showyourwork`` command and various subcommands.

All commands and subcommands are executed within an isolated conda environment
running the version of ``showyourwork`` specified by the article's ``version``
setting. This allows users running any version of ``showyourwork`` to build any
article.

"""
from .main import entry_point
