import os

import pytest


def pytest_addoption(parser):
    """Add custom pytest options."""
    parser.addoption(
        "--remote",
        action="store_true",
        dest="remote",
        default=False,
        help="enable remote tests",
    )
    parser.addoption(
        "--action-spec",
        action="store",
        default="",
        help="version spec of showyourwork to install on GH Actions",
    )
    parser.addoption(
        "--action-version",
        action="store",
        default="",
        help="version of the showyourwork-action to use in the workflow",
    )
    parser.addoption(
        "--github-org",
        action="store",
        default=None,
        help="GitHub organization for test repos (defaults to 'showyourwork')",
    )
    parser.addoption(
        "--no-org",
        action="store_true",
        default=False,
        help="Use personal GitHub account instead of an organization",
    )
    parser.addoption(
        "--require-zenodo",
        action="store_true",
        default=False,
        help="Force zenodo unit tests to run even if tokens are not set",
    )


def pytest_collection_modifyitems(config, items):
    """
    Process custom CLI flags and environment variables,
    and set markers accordingly:

    - skip remote tests when ``--remote`` is used
    - skip zenodo unit tests when SANDBOX_TOKEN is missing and
      ``--require-zenodo`` is not used
    """
    # Skip remote integration tests
    if config.getoption("--remote"):
        return
    skipper = pytest.mark.skip(reason="need --remote option to run")
    for item in items:
        if "remote" in item.keywords:
            item.add_marker(skipper)

    # Skip Zenodo unit tests
    # TODO: Make sure this does not skip integration test by mistake.
    # If it does, improve item filtering or move back to unit/conftest.py
    # Same for remote block above
    sandbox_token = os.getenv("SANDBOX_TOKEN")
    require_zenodo = config.getoption("--require-zenodo")

    if not require_zenodo and not sandbox_token:
        skipper = pytest.mark.skip(reason="SANDBOX_TOKEN environment variable not set")
        for item in items:
            if "zenodo" in item.keywords:
                item.add_marker(skipper)


def pytest_configure(config):
    """
    Set environment variable according to the following CLI options:

    - ``ACTION_SPEC`` based on ``action-spec``
    - ``ACTION_VERSION`` based on ``action-version``
    - ``SHOWYOURWORK_TEST_ORG`` based on ``--no-org`` and ``--github-org``
      (defaults to ``showyourwork``)
    """
    os.environ["ACTION_SPEC"] = str(config.getoption("--action-spec"))
    os.environ["ACTION_VERSION"] = str(config.getoption("--action-version"))

    # Set GITHUB_ORG from CLI option or environment variable
    if config.getoption("--no-org"):
        github_org = None
    else:
        github_org = config.getoption("--github-org")
        if github_org is None:
            github_org = os.getenv("SHOWYOURWORK_TEST_ORG", "showyourwork")
        if github_org.lower() in ("none", ""):
            github_org = None
    os.environ["SHOWYOURWORK_TEST_ORG"] = github_org or ""
