import argparse
import json
from pathlib import Path
from typing import Dict, Optional, Set

import graphviz

parser = argparse.ArgumentParser()
parser.add_argument("--config")
parser.add_argument("--repo-path")
parser.add_argument("--work-path")
parser.add_argument("--output")
args = parser.parse_args()

with open(args.config, "r") as f:
    config = json.load(f)

repo_path = Path(args.repo_path)
work_path = Path(args.work_path)


def relative_or_skip(input_path: str) -> Optional[str]:
    path = Path(input_path)
    if path.is_relative_to(work_path):
        return None
    try:
        return str(Path(path).relative_to(repo_path))
    except ValueError:
        return None


graph_attr = {
    "ranksep": 1,
    "nodesep": 0.65,
}
node_attr = {
    "shape": "box",
    "penwidth": 2,
    "width": 1,
}
graphviz.Digraph(
    "dag",
    node_attr=node_attr,
    graph_attr=graph_attr,
    engine="sfdp",
)

full_tree = config["_dependency_tree"]
tree: Dict[str, Set[str]] = {}
for file, parents in full_tree.items():
    path = relative_or_skip(file)
    if path is None:
        continue
    tree[str(path)] = set()
    for parent in parents:
        todo = [parent]
        while len(todo):
            query = todo.pop()
            qpath = relative_or_skip(query)
            if qpath is None:
                todo += full_tree.get(query, [])
            else:
                tree[str(path)] |= {str(qpath)}

print(tree)
