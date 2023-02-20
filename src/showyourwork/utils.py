import shutil
from pathlib import Path
from typing import Optional

from showyourwork.paths import PathLike, path_to_rule_name


def copy_file_or_directory(src: PathLike, dst: PathLike) -> None:
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    if Path(src).is_dir():
        shutil.copytree(src, dst)
    else:
        shutil.copyfile(src, dst)


def rule_name(
    *parts: str, document: Optional[PathLike] = None, plugin: Optional[str] = None
) -> str:
    # The prefix is used to fake namespaces for rules
    if plugin is None:
        prefix = "syw"
    else:
        if plugin.startswith("showyourwork.plugins."):
            plugin_name = plugin[len("showyourwork.plugins.") :]
        else:
            plugin_name = plugin
        prefix = f"sywplug_{plugin_name.replace('.', '_')}"

    # The suffix is used to support multiple documents in a single project
    suffix = ""
    if document is not None:
        suffix = f"__{path_to_rule_name(document)}"

    return f"{prefix}__{'_'.join(parts)}{suffix}"
