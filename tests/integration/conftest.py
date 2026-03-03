import os
import shutil
from pathlib import Path

import pytest


def pytest_sessionstart(session):
    # Clean the sandbox
    for folder in (Path(__file__).parents[0] / "sandbox").glob("*"):
        if folder.is_dir():
            shutil.rmtree(folder)


def pytest_addoption(parser):
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


def pytest_configure(config):
    config.addinivalue_line("markers", "remote: a test that requires remote access")
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


def pytest_collection_modifyitems(config, items):
    if config.getoption("--remote"):
        return
    skipper = pytest.mark.skip(reason="need --remote option to run")
    for item in items:
        if "remote" in item.keywords:
            item.add_marker(skipper)
