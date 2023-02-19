import argparse
import json
from pathlib import Path

import graphviz

parser = argparse.ArgumentParser()
parser.add_argument("--config")
parser.add_argument("--repo-path")
parser.add_argument("--work-path")
parser.add_argument("--thumbnails-path")
parser.add_argument("--output")
args = parser.parse_args()

with open(args.config, "r") as f:
    config = json.load(f)

repo_path = Path(args.repo_path)
work_path = Path(args.work_path)
thumbnails_path = Path(args.thumbnails_path)

tree = config["_dependency_tree_simple"]
all_files = list(sorted(set(tree.keys()) | set().union(*tree.values())))

# Set up the graphviz graph for the DAG
dag_config = config.get("dag", {})
graph_attr = dag_config.get(
    "graph_attr",
    {
        "ranksep": "1",
        "nodesep": "0.65",
    },
)
node_attr = dag_config.get(
    "node_attr",
    {
        "shape": "box",
        "penwidth": "2",
        "width": "1",
    },
)
dot = graphviz.Digraph(
    "dag",
    node_attr=node_attr,
    graph_attr=graph_attr,
    engine=dag_config.get("engine", "sfdp"),
)

with dot.subgraph(name="files") as c:
    c.attr(
        color="black",
        fillcolor="#2ca02c",
        style="filled",
        label="files",
    )
    for label in all_files:
        thumbnail = thumbnails_path / label
        thumbnail = thumbnail.with_name(thumbnail.name + ".png")
        if thumbnail.is_file():
            c.node(
                label,
                label="",
                image=str(thumbnail),
                penwidth="0",
                fixedsize="true",
                imagescale="true",
                height="1",
                width="1",
                color="black",
                style="filled",
                fillcolor="white",
            )
        else:
            c.node(
                label,
                label=Path(label).name,
                # URL=file2url(file),
                color="black",
                style="filled",
                fillcolor="white",
            )


for file in all_files:
    for dependency in tree.get(file, []):
        dot.edge(dependency, file, color="black")

# Save the graph
output_path = Path(args.output)
assert output_path.suffix == ".pdf"
dot.render(directory=str(output_path.parent), filename=output_path.with_suffix("").name)
