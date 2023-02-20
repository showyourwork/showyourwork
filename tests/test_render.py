import platform

import pytest
from showyourwork import test_util


@pytest.mark.skipif(
    platform.system() == "Windows",
    reason="imagemagick can't be installed on Windows using conda",
)
def test_render() -> None:
    test_util.run("tests/projects/render", snakemake_args=["dag.pdf"])
