import os
import re
import shutil
from collections import namedtuple
from contextlib import contextmanager
from textwrap import TextWrapper

import click

from .. import __version__, exceptions, git
from . import commands

SUBCOMMANDS = ["build", "cache", "clean", "setup", "tarball"]
OPTIONS = ["-v", "--version", "--help"]
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


# Common options:
option_spec = namedtuple("option_spec", ["args", "kwargs"])
cores_option = option_spec(
    args=["-c", "--cores"],
    kwargs=dict(
        default="1",
        help="Number of cores to use; passed to snakemake.",
    ),
)
conda_frontend_option = option_spec(
    args=["--conda-frontend"],
    kwargs=dict(
        default="conda",
        help="The conda frontend to use; passed to snakemake.",
    ),
)
project_option = option_spec(
    args=["-p", "--project"],
    kwargs=dict(
        default=".",
        help="The root folder of the showyourwork project.",
    ),
)
branch_option = option_spec(
    args=["-b", "--branch"],
    kwargs=dict(
        required=False,
        default=None,
        help="Branch whose deposit is to be deleted. Default is current branch.",
    ),
)


def main_command(f):
    return main.command(context_settings=dict(ignore_unknown_options=True))(f)


def ensure_top_level():
    """Ensures we're running commands in the top level of a git repo.

    Raises:
        exceptions.ShowyourworkException: If not in the top level of a git repo.
    """
    root = os.path.realpath(git.get_repo_root())
    here = os.path.realpath(".")
    if not root == here:
        raise exceptions.ShowyourworkException(
            "The `showyourwork` command must be called "
            "from the top level of a git repository."
            "\n"
            "`showyourwork` by itself runs a document build. "
            "Did you mean to run `showyourwork setup`?"
        )


@contextmanager
def cwd_as(path):
    """Temporarily change the current working directory to `path`."""
    old_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)


def echo(text="", **kwargs):
    """Print a message to the terminal with some custom formatting.

    Breaks long lines using ``TextWrapper`` and adds custom colors and
    indentation to code snippets and environment variables.

    This may just be the ugliest function I've ever written in my life.

    Args:
        text (str): The text to print.
        kwargs: Additional keyword arguments to pass to ``click.echo``.
    """
    try:
        terminal_size = shutil.get_terminal_size().columns
    except Exception:
        terminal_size = 80
    wrapper = TextWrapper(
        width=terminal_size,
        drop_whitespace=True,
        initial_indent="",
    )
    text = text.replace("\n", " ")
    text = re.sub(r"``(.*?)``\s*", r"<SPLIT><BR><TAB>`\1`<BR><SPLIT>", text)
    lines = [line.strip() for line in text.split("<SPLIT>") if line.strip() != ""]
    for n, original_text in enumerate(lines):
        text = original_text
        while "  " in text:
            text = text.replace("  ", " ")
        text = wrapper.fill(text)
        L, R = click.style("><", fg="blue").split("><")
        text = re.sub("`(.*?)`", L + r"\1" + R, text)
        text = text.replace("<TAB>", "    ")
        text = text.replace("<BR>", "\n")
        if n == len(lines) - 1:
            if text.endswith("\n"):
                text = text[:-1]
        click.echo(text, **kwargs)


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(
    __version__,
    "--version",
    "-v",
    package_name="showyourwork",
    message="%(version)s",
)
def main():
    """Easily build open-source, reproducible scientific articles."""
    pass


@main_command
@click.option(*cores_option.args, **cores_option.kwargs)
@click.option(*conda_frontend_option.args, **conda_frontend_option.kwargs)
@click.option(*project_option.args, **project_option.kwargs)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def build(cores, conda_frontend, project, snakemake_args):
    """Build an article in the current working directory."""
    with cwd_as(project):
        ensure_top_level()
        commands.preprocess(
            snakemake_args=snakemake_args,
            cores=cores,
            conda_frontend=conda_frontend,
        )
        commands.build(
            snakemake_args=snakemake_args,
            cores=cores,
            conda_frontend=conda_frontend,
        )


def validate_slug(context, param, slug):
    """Checks whether the repo slug was provided correctly."""

    def pause():
        if context.params.get("yes"):
            pass
        else:
            click.echo("\nPress any key to continue, or Ctrl+C to abort...")
            click.getchar()
            click.echo()

    if "/" in slug and len(slug.split("/")) == 2:
        user, repo = slug.split("/")

        if not context.params.get("quiet"):
            # Greeting
            echo(
                f"""
                Let's get you set up with a new repository. I'm going to create
                a folder called ``{repo}`` in the current working directory. If
                you haven't done this yet, please visit
                ``https://github.com/new`` at this time and create an empty
                repository called ``{slug}``
                """
            )
            pause()

            # Check Zenodo credentials
            cache = context.params.get("cache")
            if not cache:
                echo(
                    """
                    By default, showyourwork does not set up remote caching for
                    this repository. To enable remote caching of datasets to
                    Zenodo, cancel this run and re-run `showyourwork setup` with
                    the `--cache` command-line option.
                    """
                )
            else:
                echo(
                    f"""
                    You requested remote caching, so I'm
                    going to create a deposit draft on Zenodo Sandbox
                    where intermediate results will be cached.
                    Please make sure at this time that you have defined the
                    `SANDBOX_TOKEN` environment variable
                    containing your API key for the service. If you don't
                    have this key, visit
                    ``https://sandbox.zenodo.org/account/settings/applications/tokens/new``
                    to create a new personal access token with
                    `deposit:actions` and `deposit:write` scopes and store it
                    in the environment variable `SANDBOX_TOKEN`.
                    In order for caching to work
                    on GitHub Actions, you'll also have to visit
                    ``https://github.com/{slug}/settings/secrets/actions/new``
                    at this time to create a corresponding repository secret
                    with the same name containing the API key.
                    """
                )
            pause()

            # Check Overleaf credentials
            if not context.params.get("overleaf"):
                echo(
                    """
                    By default, showyourwork does not set up Overleaf
                    integration. To enable Overleaf integration with this
                    repository, cancel this run and re-run `showyourwork setup`
                    with the `--overleaf` command-line option and pass in your
                    Overleaf project ID.
                    """
                )
            else:
                echo(
                    f"""
                    You provided an Overleaf project id, so I'm going to set up
                    Overleaf integration for this repository. Please make sure
                    at this time that you have defined the `OVERLEAF_TOKEN`
                    environment variable. In order for
                    this to work on GitHub Actions, please go to
                    ``https://github.com/{slug}/settings/secrets/actions/new``
                    at this time and create `OVERLEAF_TOKEN` secret with
                    your Overleaf token.
                    """
                )
            pause()

        return slug

    else:
        raise click.BadParameter("Must have the form `user/repo`.")


@main.command()
@click.argument("slug", callback=validate_slug, metavar="<user/repo>")
@click.option(
    "-y",
    "--yes",
    is_flag=True,
    default=False,
    help="Accept all `Press any key to continue` prompts automatically.",
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    default=False,
    help="Don't prompt the user, and don't display informational output.",
)
@click.option(
    "-c",
    "--cache",
    is_flag=True,
    default=False,
    help="Set up intermediate result caching on Zenodo Sandbox. "
    "Requires a Zenodo Sandbox API token provided as the "
    "environment variable and GitHub repository secret `SANDBOX_TOKEN`.",
)
@click.option(
    "-o",
    "--overleaf",
    help="Overleaf project id to sync with (optional). Requires Overleaf "
    "credentials, provided as the environment variables and GitHub repository "
    "secrets `OVERLEAF_TOKEN`.",
    default=None,
    type=type("PROJECT_ID", (str,), {}),
)
@click.option(
    "-s",
    "--ssh",
    is_flag=True,
    help="Use ssh to authenticate with GitHub? Default is to use https.",
)
@click.option(
    "-a",
    "--action-spec",
    help="Version spec of showyourwork to use in the GitHub action. Default is the "
    "latest version on PyPI.",
    default=None,
    type=type("SPEC", (str,), {}),
)
def setup(slug, yes, quiet, cache, overleaf, ssh, action_spec):
    """
    Set up a new article repository in the current working directory.

    This command expects a single positional argument, `SLUG`, of the form
    `user/repo`, where `user` is the user's GitHub handle and `repo` is the
    name of the repository (and local directory) to create.
    """
    commands.setup(slug, cache, overleaf, ssh, action_spec)


@main_command
@click.option(*cores_option.args, **cores_option.kwargs)
@click.option(*conda_frontend_option.args, **conda_frontend_option.kwargs)
@click.option(*project_option.args, **project_option.kwargs)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Forcefully remove everything in the `src/tex/figures` and `src/data` "
    "directories.",
)
@click.option(
    "-d",
    "--deep",
    is_flag=True,
    help="Forcefully remove the `.snakemake` and `~/.showyourwork` directories.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def clean(cores, conda_frontend, project, force, deep, snakemake_args):
    """Clean the article build in the current working directory."""
    with cwd_as(project):
        ensure_top_level()
        commands.clean(
            force,
            deep,
            snakemake_args=snakemake_args,
            cores=cores,
            conda_frontend=conda_frontend,
        )


@main_command
@click.option(*cores_option.args, **cores_option.kwargs)
@click.option(*conda_frontend_option.args, **conda_frontend_option.kwargs)
@click.option(*project_option.args, **project_option.kwargs)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def tarball(cores, conda_frontend, project, snakemake_args):
    """Generate a tarball of the build in the current working directory."""
    with cwd_as(project):
        ensure_top_level()
        commands.preprocess(
            snakemake_args=snakemake_args,
            cores=cores,
            conda_frontend=conda_frontend,
        )
        commands.tarball(
            snakemake_args=snakemake_args,
            cores=cores,
            conda_frontend=conda_frontend,
        )


@main.group()
@click.pass_context
def cache(ctx):
    """Caching-related operations."""
    pass


@cache.command()
@click.pass_context
@click.option(*branch_option.args, **branch_option.kwargs)
def create(ctx, branch):
    """
    Create a Zenodo Sandbox deposit draft for the given branch.

    Requires a Zenodo Sandbox API token provided as the
    environment variable and Github repository secret `SANDBOX_TOKEN`.
    """
    ensure_top_level()
    commands.zenodo_create(branch)


@cache.command()
@click.pass_context
@click.option(*branch_option.args, **branch_option.kwargs)
def delete(ctx, branch):
    """
    Delete the latest draft of a Zenodo Sandbox deposit.

    Requires a Zenodo Sandbox API token provided as the
    environment variable and Github repository secret `SANDBOX_TOKEN`.
    """
    ensure_top_level()
    commands.zenodo_delete(branch)


@cache.command()
@click.pass_context
@click.option(*branch_option.args, **branch_option.kwargs)
def freeze(ctx, branch):
    """
    Publishes the current Zenodo Sandbox deposit draft for the given branch
    to Zenodo Sandbox.

    Useful for preserving the current cache state even when changes are made
    to the rule inputs. Note that Zenodo Sandbox storage is temporary, so if
    you'd like to publish your cache on Zenodo, you should call ``publish``
    instead.

    Requires a Zenodo Sandbox API token provided as the
    environment variable and Github repository secret `SANDBOX_TOKEN`.
    """
    ensure_top_level()
    commands.zenodo_freeze(branch)


@cache.command()
@click.pass_context
@click.option(*branch_option.args, **branch_option.kwargs)
def publish(ctx, branch):
    """
    Publishes the current Zenodo Sandbox deposit draft for the given branch
    to Zenodo for permanent, public storage. This action assigns a permanent,
    static DOI to the deposit, which cannot be deleted or undone.

    Requires a Zenodo API token provided as the environment variable
    `ZENODO_TOKEN`.
    """
    ensure_top_level()
    commands.zenodo_publish(branch)


@cache.command(hidden=True)  # used internally
@click.pass_context
def restore(ctx):
    ensure_top_level()
    commands.cache_restore()


@cache.command(hidden=True)  # used internally
@click.pass_context
def update(ctx):
    ensure_top_level()
    commands.cache_update()
