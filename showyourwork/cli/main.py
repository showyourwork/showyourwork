from . import commands
from .. import git
from .. import exceptions
import os
import click


@click.group(invoke_without_command=True)
@click.pass_context
def entry_point(context):
    """Easily build open-source, reproducible scientific articles."""
    if context.invoked_subcommand is None:
        # Default command is `build`
        context.invoke(build)
    elif context.invoked_subcommand == "setup":
        pass
    else:
        # Ensure we're running the command in the top level of a git repo
        root = os.path.realpath(git.get_repo_root())
        here = os.path.realpath(".")
        if not root == here:
            raise exceptions.ShowyourworkException(
                "Must run `showyourwork` in the top level of a git repo."
            )


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


@entry_point.command(hidden=True)
@click.option("--restore", is_flag=True)
@click.option("--update", is_flag=True)
def cache(restore, update):
    """Update or restore the cache on GitHub Actions."""
    if not restore != update:
        raise exceptions.ShowyourworkException(
            "Must provide either `--restore` or `--update`."
        )
    if restore:
        commands.cache_restore()
    else:
        commands.cache_update()