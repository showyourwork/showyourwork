import json

from showyourwork.paths import PathLike


def write_manuscript_dependencies(input_file: PathLike, output_file: PathLike) -> None:
    with open(input_file, "r") as f:
        dependencies = json.load(f)
    files = list(dependencies.get("unlabled", [])) + list(dependencies.get("files", []))
    for figure in dependencies.get("figures", {}).values():
        files.extend(figure)
    with open(output_file, "w") as f:
        json.dump(files, f, indent=2)
