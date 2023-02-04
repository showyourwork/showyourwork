from importlib.resources import as_file, files
from pathlib import Path
from typing import Any, Dict, Union

PathLike = Union[str, Path]


def package_data(module: str, *file: str) -> Path:
    with as_file(files(module).joinpath(*file)) as f:
        path = Path(f)
    path = path.resolve()
    if not path.exists():
        raise FileNotFoundError(
            f"No file exists at the path {file} for module {module}. This might have "
            "something to do with how you installed showyourwork or a plugin. Open an "
            "issue on the relevant repository."
        )
    return path


class PathMeta:
    root: Path

    def __truediv__(self, name: PathLike) -> Path:
        return self.subdir(name)

    def subdir(self, name: PathLike, *others: PathLike) -> Path:
        path = self.root / name
        for other in others:
            path = path / other
        path.mkdir(parents=True, exist_ok=True)
        return path


class repo(PathMeta):
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.root = find_project_root()

    @property
    def manuscript(self) -> Path:
        return self.root / self.config.get("manuscript", "ms.tex")


class work(PathMeta):
    def __init__(self, config: Dict[str, Any]):
        working_directory = config.get("working-directory", None)
        if working_directory is None:
            working_directory = find_project_root()
        else:
            working_directory = Path(working_directory)
        working_directory = working_directory / ".showyourwork"
        working_directory.mkdir(parents=True, exist_ok=True)
        self.root = working_directory

    def plugin(self, name: str, *others: PathLike) -> Path:
        return self.subdir("plugins", name, *others)


def find_project_root(*input_paths: PathLike) -> Path:
    """Find the root of the project, defined as the first parent directory that
    contains a .git directory or showyourwork.yml file. Based on the implementation from
    psf/black.
    """
    if input_paths:
        srcs = list(input_paths)
    else:
        srcs = [Path.cwd().resolve()]
    path_srcs = [Path(Path.cwd(), src).resolve() for src in srcs]

    src_parents = [
        list(path.parents) + ([path] if path.is_dir() else []) for path in path_srcs
    ]
    common_base = max(
        set.intersection(*(set(parents) for parents in src_parents)),
        key=lambda path: path.parts,
    )

    for directory in (common_base, *common_base.parents):
        if (directory / ".git").exists():
            return directory
        if (directory / "showyourwork.yml").is_file():
            return directory
        if (directory / "showyourwork.yaml").is_file():
            return directory

    return directory
