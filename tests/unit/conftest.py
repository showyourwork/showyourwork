import os

import pytest


def pytest_addoption(parser):
    """Add custom pytest options."""
    parser.addoption(
        "--require-zenodo",
        action="store_true",
        default=False,
        help="Force zenodo tests to run even if tokens are not set",
    )


def pytest_configure(config):
    """Register custom pytest markers."""
    config.addinivalue_line(
        "markers",
        "zenodo: a test that requires Zenodo/Zenodo Sandbox API access",
    )


def pytest_collection_modifyitems(config, items):
    """
    Skip zenodo tests if SANDBOX_TOKEN is not set (unless --require-zenodo is passed).
    """
    sandbox_token = os.getenv("SANDBOX_TOKEN")
    zenodo_token = os.getenv("ZENODO_TOKEN")
    require_zenodo = config.getoption("--require-zenodo")

    if not require_zenodo and (not sandbox_token or not zenodo_token):
        skipper = pytest.mark.skip(
            reason="SANDBOX_TOKEN or ZENODO_TOKEN environment variable not set"
        )
        for item in items:
            if "zenodo" in item.keywords:
                item.add_marker(skipper)
