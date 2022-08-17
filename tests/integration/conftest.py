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


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "remote: a test that requires remote access"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--remote"):
        return
    skipper = pytest.mark.skip(reason="need --remote option to run")
    for item in items:
        if "remote" in item.keywords:
            item.add_marker(skipper)
