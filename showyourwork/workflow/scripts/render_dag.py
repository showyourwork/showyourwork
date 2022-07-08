"""
Generates a directed acyclic graph (DAG) of the build process.

"""
from showyourwork import paths, exceptions
from showyourwork.subproc import get_stdout
from showyourwork.zenodo import get_dataset_dois
from pathlib import Path
from jinja2 import Environment, BaseLoader
import graphviz
import json

# Convert to PNG & get aspect ratio
def convert_to_png(file):
    """
    Convert a file to a PNG thumbnail and return its aspect ratio.

    On error, returns None.
    """

    def get_aspect(code, stdout, stderr):
        try:
            if code != 0:
                raise exceptions.CalledProcessError(stderr)
            aspect = float(
                get_stdout(
                    f'convert {file} -format "%[fx:w/h]" info:', shell=True,
                )
            )
        except exceptions.CalledProcessError:
            return None
        else:
            return aspect

    aspect = get_stdout(
        f"convert -resize 500x500 {file} {file}.png",
        shell=True,
        callback=get_aspect,
    )

    return aspect


if __name__ == "__main__":

    # Snakemake config (available automagically)
    config = snakemake.config  # type:ignore

    # Colors
    colors = dict(
        zenodo="#1f77b4",
        script="#1FB45B",
        tex="black",
        dataset="#AC5678",
        other="#636363",
        figures="black",
        edge="black",
    )

    # Instantiate the graph
    dot = graphviz.Digraph(
        "dag",
        node_attr={"shape": "box", "penwidth": "2", "width": "1"},
        engine="sfdp",
    )

    # Collect all dependency information
    dependencies = config["dag_dependencies"]

    # Ignore temporary & showyourwork files
    ignore = config["tex_files_out"] + [
        config["stylesheet"],
        config["stylesheet_meta_file"],
        str(paths.user().flags / "SYW__DAG"),
        str(paths.user().compile),
        "dag.pdf",
        "showyourwork.yml",
    ]

    # Assemble all files
    files = []
    for file in config["dag_dependencies"]:
        files.append(file)
        files.extend(config["dag_dependencies"][file])
    files = list(set(files))
    files = [file for file in files if file not in ignore]

    # Group into different types
    datasets = []
    scripts = []
    texfiles = []
    figures = []
    others = []
    for file in files:
        if Path(file).absolute().is_relative_to(paths.user().data):
            datasets.append(file)
        elif Path(file).absolute().is_relative_to(paths.user().scripts):
            scripts.append(file)
        elif Path(file).absolute().is_relative_to(paths.user().figures):
            figures.append(file)
        elif Path(file).absolute().is_relative_to(paths.user().tex):
            texfiles.append(file)
        elif file == config["ms_pdf"]:
            pass
        else:
            others.append(file)

    # Zenodo nodes
    with dot.subgraph(name="cluster_zenodo") as c:
        c.attr(color=colors["zenodo"])
        # Static
        for url in get_dataset_dois(datasets, config["datasets"]):
            doi = url.split("record/")[-1]
            c.node(doi, color=colors["zenodo"])
        # Cache
        branch = config["git_branch"]
        cache = (
            config["cache"][branch]["zenodo"]
            or config["cache"][branch]["sandbox"]
        )
        if cache:
            c.node(cache, color=colors["zenodo"])

    # Dataset nodes
    with dot.subgraph(name="cluster_data") as c:
        c.attr(color=colors["dataset"])
        for file in datasets:
            label = file.removeprefix("src/data/")
            c.node(file, label=label, color=colors["dataset"])
            for url in get_dataset_dois([file], config["datasets"]):
                c.edge(doi, file)

    # Script nodes
    with dot.subgraph(name="cluster_scripts") as c:
        c.attr(color=colors["script"])
        for file in scripts:
            label = file.removeprefix("src/scripts/")
            c.node(file, label=label, color=colors["script"])

    # Tex nodes
    with dot.subgraph(name="cluster_tex") as c:
        c.attr(color=colors["tex"])
        for file in texfiles:
            label = file.removeprefix("src/tex/")
            c.node(file, label=label, color=colors["tex"])

    # Figure nodes
    with dot.subgraph(name="cluster_figures") as c:
        c.attr(color=colors["figures"])
        for file in figures:
            if aspect := convert_to_png(file):
                # Add PNG thumbnail
                c.node(
                    file,
                    label="",
                    image=f"{file}.png",
                    penwidth="0",
                    fixedsize="true",
                    imagescale="false",
                    width=f"{aspect}",
                    height="1",
                )
            else:
                # Couldn't convert to PNG; display text instead
                label = file.removeprefix("src/data/")
                c.node(file, label=label, color=colors["figures"])

    # Other nodes
    with dot.subgraph(name="cluster_others") as c:
        c.attr(color=colors["other"])
        for file in others:
            c.node(file, color=colors["other"])

    # Article node
    dot.node(
        config["ms_pdf"],
        label="",
        image=str(
            paths.showyourwork().resources / "img" / "article-thumb.png"
        ),
        penwidth="0",
        fixedsize="true",
        imagescale="false",
        width="1.5",
        height="1.5",
    )

    # Add the edges
    for file in files:
        for dependency in dependencies.get(file, []):
            if dependency not in ignore:
                dot.edge(dependency, file, color=colors["edge"])

                # Cache
                if dependency in config["cached_deps"]:
                    dot.edge(
                        cache, dependency, color=colors["edge"], style="dashed"
                    )
                if file in config["cached_deps"]:
                    dot.edge(
                        dependency, cache, color=colors["edge"], style="dashed"
                    )

    # Save the graph
    dot.save(directory=paths.user().repo)

    # Render the PDF
    get_stdout("dot -Tpdf dag.gv > dag.pdf", shell=True, cwd=paths.user().repo)

    # Remove temporary files
    for figure in figures:
        pngfile = Path(f"{figure}.png")
        if pngfile.exists():
            pngfile.unlink()
    gvfile = paths.user().repo / "dag.gv"
    if gvfile.exists():
        gvfile.unlink()
