from pathlib import Path
from typing import List


def snakefiles() -> List[Path]:
    from showyourwork.paths import package_data

    return [package_data("showyourwork.plugins.tex", "workflow", "Snakefile")]
