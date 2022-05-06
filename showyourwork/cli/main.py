from . import commands
from .. import git, exceptions, zenodo, __version__
from textwrap import TextWrapper
import os
import shutil
import click
import re


def ensure_top_level():
    """
    Ensures we're running commands in the top level of a git repo.

    """
    root = os.path.realpath(git.get_repo_root())
    here = os.path.realpath(".")
    if not root == here:
        raise exceptions.ShowyourworkException(
            "The `showyourwork` command must be called "
            "from the top level of a git repo."
        )


def echo(text="", **kwargs):
    """
    Print a message to screen with some custom formatting.

    This may be the ugliest function I've ever written in my life.

    """
    try:
        terminal_size = shutil.get_terminal_size().columns
    except:
        terminal_size = 80
    wrapper = TextWrapper(
        width=terminal_size,
        drop_whitespace=True,
        initial_indent="",
    )
    text = text.replace("\n", " ")
    text = re.sub("``(.*?)``\s*", r"<SPLIT><BR><TAB>`\1`<BR><SPLIT>", text)
    lines = [
        line.strip() for line in text.split("<SPLIT>") if line.strip() != ""
    ]
    for n, text in enumerate(lines):
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


@click.group(invoke_without_command=True)
@click.option(
    "-v",
    "--version",
    is_flag=True,
    help="Show the program version and exit.",
)
@click.pass_context
def entry_point(context, version):
    """Easily build open-source, reproducible scientific articles."""
    # Parse
    if version:
        print(__version__)
    elif context.invoked_subcommand is None:
        # Default command is `build`
        context.invoke(build)


@entry_point.command(
    context_settings=dict(
        ignore_unknown_options=True,
    )
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def build(snakemake_args):
    """Build an article in the current working directory."""
    ensure_top_level()
    commands.preprocess()
    commands.build(snakemake_args)


def validate_slug(context, param, slug):
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
            sandbox = context.params.get("sandbox")
            if not cache:
                echo(
                    """
                    You didn't request remote caching (via the
                    `--cache` command-line option), 
                    so I'm not going to set up remote caching for
                    this repository.
                    """
                )
            else:
                if sandbox:
                    service = zenodo.services["sandbox"]
                else:
                    service = zenodo.services["zenodo"]
                echo(
                    f"""
                    You requested remote caching on {service["name"]}, so I'm
                    going to create a deposit draft where intermediate results
                    will be cached. Please make sure at this time that you have
                    defined the `{service["token_name"]}` environment variable
                    containing your API key for {service["name"]}. If you don't
                    have one, visit
                    ``https://{service['url']}/account/settings/applications/tokens/new``
                    to create a new personal access token with
                    `deposit:actions` and `deposit:write` scopes and store it
                    in the environment variable `{service["token_name"]}`. In
                    order for this to work on GitHub Actions, you'll also have
                    to visit
                    ``https://github.com/{slug}/settings/secrets/actions/new``
                    at this time to create a `{service["token_name"]}` secret
                    with your API access token.
                    """
                )
            pause()

            # Check Overleaf credentials
            if not context.params.get("overleaf"):
                echo(
                    f"""
                    You didn't provide an Overleaf project id (via the 
                    `--overleaf` command-line option), so I'm not going to set 
                    up Overleaf integration for this repository.
                    """
                )
            else:
                echo(
                    f"""
                    You provided an Overleaf project id, so I'm going to set up 
                    Overleaf integration for this repository. Please make sure 
                    at this time that you have defined the `OVERLEAF_EMAIL` and 
                    `OVERLEAF_PASSWORD` environment variables. In order for 
                    this to work on GitHub Actions, please go to
                    ``https://github.com/{slug}/settings/secrets/actions/new``
                    at this time and create `OVERLEAF_EMAIL` and 
                    `OVERLEAF_PASSWORD` secrets with your Overleaf credentials.
                    """
                )
            pause()

        return slug

    else:

        raise click.BadParameter("Must have the form `user/repo`.")


@entry_point.command()
@click.argument("slug", callback=validate_slug)
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
    help="Set up intermediate result caching on Zenodo. Requires a Zenodo "
    "API token provided as the environment variable and Github repository "
    "secret ZENODO_TOKEN.",
)
@click.option(
    "-x",
    "--sandbox",
    is_flag=True,
    default=False,
    help="Use Zenodo Sandbox to cache results (if `--cache` was provided). "
    "Requires a Zenodo Sandbox API token provided as the environment variable "
    "and Github repository secret SANDBOX_TOKEN.",
)
@click.option(
    "-o",
    "--overleaf",
    help="Overleaf project id to sync with (optional). Requires Overleaf "
    "credentials, provided as the environment variables and GitHub repository "
    "secrets OVERLEAF_EMAIL and OVERLEAF_PASSWORD.",
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
    "-v",
    "--version",
    help="Version of showyourwork to use. Default is the current version.",
    default=None,
    type=type("VERSION", (str,), {}),
)
def setup(
    slug,
    yes,
    quiet,
    cache,
    sandbox,
    overleaf,
    ssh,
    version,
):
    """
    Set up a new article repository in the current working directory.

    This command expects a single positional argument, `SLUG`, of the form
    `user/repo`, where `user` is the user's GitHub handle and `repo` is the
    name of the repository (and local directory) to create.
    """
    if cache:
        if sandbox:
            cache_service = "sandbox"
        else:
            cache_service = "zenodo"
    else:
        if sandbox:
            raise exceptions.ShowyourworkException(
                "Must provide the `--cache` option if `--sandbox` is provided."
            )
        else:
            cache_service = None
    commands.setup(slug, cache_service, overleaf, ssh, version)


@entry_point.command()
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Forcefully remove everything in the `src/tex/figures` and `src/data` directories.",
)
@click.option(
    "-d",
    "--deep",
    is_flag=True,
    help="Forcefully remove the `.snakemake` and `~/.showyourwork` directories.",
)
def clean(force, deep):
    """Clean the article build in the current working directory."""
    ensure_top_level()
    commands.clean(force, deep)


@entry_point.command()
def tarball():
    """Generate a tarball of the build in the current working directory."""
    ensure_top_level()
    commands.preprocess()
    commands.tarball()


@entry_point.group()
@click.pass_context
def cache(ctx):
    pass


@cache.command()
@click.pass_context
@click.option(
    "-b",
    "--branch",
    required=False,
    default=None,
    help="Which branch to create the deposit for. Default is current branch.",
)
@click.option(
    "-x",
    "--sandbox",
    is_flag=True,
    help="Use Zenodo Sandbox for caching.",
)
def create(ctx, branch, sandbox):
    """
    Create a Zenodo deposit draft for the given branch.

    Requires a Zenodo or Zenodo Sandbox API token provided as the
    environment variable and Github repository secret ZENODO_TOKEN
    or SANDBOX_TOKEN.
    """
    ensure_top_level()
    if sandbox:
        service = "sandbox"
    else:
        service = "zenodo"
    commands.zenodo_create(branch, service)


@cache.command()
@click.pass_context
@click.option(
    "-b",
    "--branch",
    required=False,
    default=None,
    help="Branch whose deposit is to be deleted. Default is current branch.",
)
def delete(ctx, branch):
    """
    Delete the latest draft of a Zenodo deposit.

    Requires a Zenodo or Zenodo Sandbox API token provided as the
    environment variable and Github repository secret ZENODO_TOKEN
    or SANDBOX_TOKEN.
    """
    ensure_top_level()
    commands.zenodo_delete(branch)


@cache.command()
@click.pass_context
@click.option(
    "-b",
    "--branch",
    required=False,
    default=None,
    help="Branch whose deposit is to be published. Default is current branch.",
)
def publish(ctx, branch):
    """
    Publish the latest draft of a Zenodo deposit.

    Requires a Zenodo or Zenodo Sandbox API token provided as the
    environment variable and Github repository secret ZENODO_TOKEN
    or SANDBOX_TOKEN.
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