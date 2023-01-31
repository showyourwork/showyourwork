"""
Command line interface for ``showyourwork``.

This module provides the command line interface for ``showyourwork``,
consisting of the ``showyourwork`` command and various subcommands.

All commands and subcommands are executed within an isolated conda environment
running the version of ``showyourwork`` specified by the article's ``version``
setting. This allows users running any version of ``showyourwork`` to build any
article.

"""
import sys
import warnings

from .. import __version__
from .main import DEFAULT_SUBCOMMAND, OPTIONS, SUBCOMMANDS, main


def entry_point():
    """
    Modify the call from `showyourwork ...` to `showyourwork build...` if the
    user didn't explicitly provide a valid subcommand or option.

    This behavior will be deprecated in a future version of `showyourwork`,
    because it can lead to some unexpected results when invalid arguments are
    provided.

    Click can in principle handle this (the `@click.group` decorator accepts an
    `invoke_without_command` option) but if we do that I don't _think_ we can
    simultaneously set `ignore_unknown_options`. We need to be able to forward
    unknown options directly to `snakemake`.

    This hack allows users to call, e.g.,

        showyourwork --rerun-incomplete

    instead of

        showyourwork build --rerun-incomplete

    (the invocation they had to use previously). We allow this by injecting the
    `build` subcommand into `sys.argv` prior to click taking control.

    """
    if len(sys.argv) == 1 or sys.argv[1] not in SUBCOMMANDS + OPTIONS:
        warnings.warn(
            "The use of the `showyourwork` command line without an explicit "
            "subcommand is deprecated. Use `showyourwork build` for "
            "consistent results.",
            FutureWarning,
        )
        sys.argv.insert(1, DEFAULT_SUBCOMMAND)

    if sys.argv[1] in ["-v", "--version"]:
        print(__version__)

    else:
        main()
