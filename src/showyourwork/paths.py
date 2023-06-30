import hashlib
from functools import cached_property, lru_cache
from importlib.resources import as_file, files
from pathlib import Path
from typing import Any, Dict, Union

PathLike = Union[str, Path]


@lru_cache
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

    raise RuntimeError(
        "Could not find project root; are you sure that you're calling showyourwork "
        "from within a showyourwork project?"
    )


def package_data(module: str, *file: str, check: bool = True) -> Path:
    with as_file(files(module).joinpath(*file)) as f:
        path = Path(f)
    path = path.resolve()
    if check and not path.exists():
        raise FileNotFoundError(
            f"No file exists at the path {'/'.join(file)} for module {module}. This "
            "might have something to do with how you installed showyourwork or a "
            "plugin. Open an issue on the relevant repository."
        )
    return path


def path_to_identifier(path: PathLike) -> str:
    path_hash = hashlib.md5(str(path).encode()).hexdigest()
    return f"{path_hash}_{Path(path).name}"


def path_to_rule_name(path: PathLike) -> str:
    return path_to_identifier(path).replace(".", "_")


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


class work(PathMeta):
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        working_directory = config.get("working_directory", None)
        if working_directory is None:
            working_directory = find_project_root() / ".showyourwork"
        else:
            working_directory = Path(working_directory)
        working_directory.mkdir(parents=True, exist_ok=True)
        self.root = working_directory

    def plugin(self, name: str, *others: PathLike) -> Path:
        return self.subdir("plugins", name, *others)

    def flag(self, name: str) -> Path:
        return self.subdir("flags") / name

    def dependencies_for(self, name: str) -> Path:
        return self.root / f"{name}.dependencies.json"

    @cached_property
    def logs(self) -> Path:
        return self.subdir("logs")

    @property
    def dependencies(self) -> Path:
        return self.root / "dependencies.json"

    @cached_property
    def build(self) -> Path:
        return self.subdir("build")
