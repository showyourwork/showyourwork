import re

from showyourwork import paths


def add_to_preamble(
    preamble_file: paths.PathLike,
    input_file: paths.PathLike,
    output_file: paths.PathLike,
) -> None:
    with open(preamble_file, "r") as f:
        preamble = f.read()
    with open(input_file, "r") as f:
        document = f.read()
    with open(output_file, "w") as f:
        f.write(_add_to_preamble(preamble, document))


def _add_to_preamble(preamble: str, document: str) -> str:
    def replace(match: re.Match[str]) -> str:
        result = match.group(0) + "\n\n" + preamble + "\n"
        return result

    return re.sub(r"(\\documentclass(?:\[.*\])?\{.*\})", replace, document)
