import sys
import warnings

import pytest


def test_requests_import_no_dependency_warning():
    """Test that importing requests doesn't produce RequestsDependencyWarning."""
    # Remove requests from sys.modules to force reimport
    modules_to_remove = [key for key in sys.modules if key.startswith("requests")]
    for module in modules_to_remove:
        del sys.modules[module]

    # Use warnings.catch_warnings to capture all warnings
    with warnings.catch_warnings(record=True) as warning_list:
        warnings.simplefilter("always")
        import requests  # noqa: F401

    # Check if any RequestsDependencyWarning was raised
    for warning in warning_list:
        if "RequestsDependencyWarning" in str(warning.category):
            pytest.fail(f"RequestsDependencyWarning detected: {warning.message}")
        if "RequestsDependencyWarning" in str(warning.message):
            pytest.fail(f"RequestsDependencyWarning detected: {warning.message}")
