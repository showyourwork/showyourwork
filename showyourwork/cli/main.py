from . import commands
import click


@click.group(invoke_without_command=True)
@click.pass_context
def entry_point(context):
    """Easily build open-source, reproducible scientific articles."""
    # Default to build
    if context.invoked_subcommand is None:
        context.invoke(build)


@entry_point.command()
def build():
    """Build an article in the current working directory."""
    commands.preprocess()
    commands.build()


@entry_point.command()
def setup():
    """Set up a new article repository in the current working directory."""
    commands.setup()


@entry_point.command()
def clean():
    """Clean the article build in the current working directory."""
    commands.clean()


@entry_point.command()
def tarball():
    """Generate a tarball of the build in the current working directory."""
    commands.preprocess()
    commands.tarball()