import platform

import pytest
from showyourwork.testing import run_showyourwork


@pytest.mark.skipif(
    platform.system() == "Windows",
    reason="imagemagick can't be installed on Windows using conda",
)
def test_render() -> None:
    run_showyourwork("tests/projects/render", "dag.pdf")
