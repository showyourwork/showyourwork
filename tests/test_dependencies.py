import json
from pathlib import Path
from typing import List

from showyourwork.dependencies import simplify_dependency_tree
from showyourwork.testing import run_showyourwork


def test_dependency_tree() -> None:
    with run_showyourwork(
        "tests/projects/dependency_tree", "syw__save_dependencies", show_diff=True
    ) as d:
        with open(d / ".showyourwork" / "dependency_tree.json", "r") as f:
            data = json.load(f)
        deps: List[str] = []
        for c in "abcdefgh":
            assert data[c] == deps
            deps.append(c)


def test_explicit_dependencies() -> None:
    run_showyourwork("tests/projects/explicit_dependencies")


work_path = Path("/work")
repo_path = Path("/repo")


def work(*x: str) -> str:
    return str(Path(work_path, *x))


def repo(*x: str) -> str:
    return str(Path(repo_path, *x))


def test_simplify() -> None:
    input_tree = {
        repo("a"): [work("b"), repo("c")],
        work("b"): [repo("c"), work("d")],
        work("d"): [repo("d"), work("e")],
        work("e"): [repo("c")],
    }
    output_tree = simplify_dependency_tree(input_tree, repo_path, work_path)
    assert output_tree == {"a": ["c", "d"]}


def test_simplify_with_directory() -> None:
    input_tree = {
        repo("a"): [work("b", "c"), repo("c")],
        work("b"): [repo("d")],
    }
    output_tree = simplify_dependency_tree(input_tree, repo_path, work_path)
    assert output_tree == {"a": ["c", "d"]}
