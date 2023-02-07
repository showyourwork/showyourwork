import json
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Generator

from showyourwork import paths, test_util


@contextmanager
def dump_dependencies(dependencies: Any) -> Generator[Path, None, None]:
    with test_util.temporary_project() as d:
        p = paths.work({"working_directory": d})
        with open(p.dependencies, "w") as f:
            json.dump(dependencies, f)
        test_util.run_snakemake(
            str(
                paths.package_data(
                    "showyourwork", "workflow", "rules", "dependencies.smk"
                )
            ),
            ["syw__dump_dependencies", "--config", f"working_directory={d}"],
            cwd=d,
        )
        deps = Path(d) / "dependency_tree.json"
        assert deps.is_file()
        yield deps


def test_empty() -> None:
    with dump_dependencies({}) as f:
        data = json.load(open(f, "r"))
        data = {Path(k).name: v for k, v in data.items()}
        assert "dependency_tree.json" in data


def test_basic() -> None:
    # This one is a little hacky. Since we don't have rules to generate these
    # files, our call to snakemake will fail. Then we're parsing the errors
    # reported by snakemake to identify the missing files, which should be the
    # dependencies we're asking for.
    #
    # TODO(dfm): It would be nice to make this a bit less brittle since it
    # currently depends on the very specific format of the error message, which
    # might break with future versions of snakemake.
    try:
        with dump_dependencies(
            {
                "unlabled": ["a", "b"],
                "files": ["c", "d"],
                "figures": {"x": ["e", "f"], "y": ["g", "h"]},
            }
        ):
            pass
    except RuntimeError as e:
        txt = str(e)
        assert "MissingInputException" in txt
        idx = txt.find("affected files:")
        deps = [f.strip() for f in txt[idx:].split("\n\n")[0].splitlines()[1:]]
        for d in ["a", "b", "c", "d", "e", "f", "g", "h"]:
            assert d in deps
