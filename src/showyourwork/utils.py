import shutil
from pathlib import Path

from showyourwork.paths import PathLike


def copy_file_or_directory(src: PathLike, dst: PathLike) -> None:
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    if Path(src).is_dir():
        shutil.copytree(src, dst)
    else:
        shutil.copyfile(src, dst)
