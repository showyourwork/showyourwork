from pathlib import Path
from typing import List

from showyourwork.plugins.staging.config import configure as configure
from showyourwork.plugins.staging.stages import NoOpStage as NoOpStage
from showyourwork.plugins.staging.stages import Stage as Stage
from showyourwork.plugins.staging.zenodo import ZenodoStage as ZenodoStage


def snakefiles() -> List[Path]:
    from showyourwork.paths import package_data

    return [package_data("showyourwork.plugins.staging", "workflow", "Snakefile")]
