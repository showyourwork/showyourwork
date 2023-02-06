from pathlib import Path
from typing import List

from showyourwork.plugins.tex.tex import add_to_preamble as add_to_preamble


def snakefiles() -> List[Path]:
    from showyourwork.paths import package_data

    return [package_data("showyourwork.plugins.tex", "workflow", "Snakefile")]
