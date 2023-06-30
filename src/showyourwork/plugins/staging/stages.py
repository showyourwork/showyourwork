from abc import ABC, abstractmethod
from collections import OrderedDict
from pathlib import Path
from typing import List, Optional

from showyourwork.paths import PathLike, package_data, path_to_identifier
from showyourwork.plugins.staging.config import _CONFIG

STAGES: OrderedDict[str, "Stage"] = OrderedDict()


class Stage(ABC):
    def __init__(
        self, name: str, restore: bool, working_directory: Optional[PathLike] = None
    ):
        self.name = name
        self.files: OrderedDict[str, PathLike] = OrderedDict()
        self.restore = restore
        self._working_directory = working_directory
        if self.name in STAGES:
            raise ValueError(f"A stage called {self.name} already exists")
        STAGES[self.name] = self

    @property
    def working_directory(self) -> Path:
        if self._working_directory is None:
            return Path(_CONFIG.get("working_directory", "staging"))
        return Path(self._working_directory)

    @property
    def upload_flag_file(self) -> Path:
        return self.working_directory / f"{self.name}.upload"

    def __call__(self, *files: PathLike) -> List[PathLike]:
        return self.staged(*files)

    def staged(self, *files: PathLike) -> List[PathLike]:
        # Add this list of files to the database of files for this stage
        for file in files:
            identifier = path_to_identifier(file)
            if identifier in self.files:
                raise RuntimeError(
                    f"Duplicate file detected in stage {self.name}: {file}"
                )
            self.files[identifier] = file

        # If we're restoring this stage, we short-circuit this rule so that it won't
        # ever be run to generate the staged files
        if self.restore:
            return []
        else:
            return list(files)

    @abstractmethod
    def snakefile(self) -> Path:
        ...


class NoOpStage(Stage):
    @property
    def directory(self) -> Path:
        return self.working_directory / self.name

    def snakefile(self) -> Path:
        return package_data("workflow", "rules", "noop.smk")
