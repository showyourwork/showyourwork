from .setup import _setup
from .build import _build
import click


@click.group()
def entry_point():
    """Command-line entry point."""
    pass


@entry_point.command()
def build():
    """Build the article."""
    _build()