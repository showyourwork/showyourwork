import json
from pathlib import Path

from showyourwork import test_util
from showyourwork.dependencies import simplify_dependency_tree


def test_dependency_tree() -> None:
    with test_util.run_context(
        "tests/projects/dependency_tree", snakemake_args=["syw__dump_dependencies"]
    ) as d:
        with open(d / ".showyourwork" / "dependency_tree.json", "r") as f:
            data = json.load(f)
        for c in "abcdefgh":
            assert data[c] == []


def test_explicit_dependencies() -> None:
    test_util.run("tests/projects/explicit_dependencies")


work_path = Path("work")
repo_path = Path("repo")


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
