import inspect
from typing import Any, Dict, List

from showyourwork.paths import PathLike

STAGES: Dict[str, List[PathLike]] = {}

DEFAULT = "default"
ARTIFACTS = "artifacts"


def get_stages_to_restore(config: Dict[str, Any]) -> List[str]:
    restore = config.get("restore", [])
    if isinstance(restore, str):
        restore = [restore]
    return restore


def staged(*files: PathLike, stage: str = DEFAULT) -> List[PathLike]:
    # Find the config in the callers frame
    frame = inspect.currentframe()
    assert frame is not None and frame.f_back is not None
    config: Dict[str, Any] = frame.f_back.f_globals["config"]

    # Extract the list of stages we want to restore from the config
    restore = get_stages_to_restore(config)

    # Add this list of files to the global database of files for this stage
    if stage not in STAGES:
        STAGES[stage] = []
    STAGES[stage].extend(files)

    # If we're restoring this stage, we short-circuit this rule so that it won't
    # ever be run to generate the staged files
    if stage in restore:
        return []
    else:
        return list(files)


def optionally_require_zenodo(config: Dict[str, Any]) -> List[PathLike]:
    return []


def snapshot_stage(config: Dict[str, Any], stage: str) -> None:
    pass


def restore_stage(config: Dict[str, Any], stage: str) -> None:
    pass
