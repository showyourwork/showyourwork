"""
Generates a directed acyclic graph (DAG) of the build process.

"""
from pathlib import Path

from jinja2 import BaseLoader, Environment

from showyourwork import exceptions, paths
from showyourwork.subproc import get_stdout
from showyourwork.zenodo import get_dataset_dois

try:
    import graphviz
except:
    graphviz = None

# Convert to PNG & get aspect ratio
def convert_to_png(file):
    """
    Convert a file to a PNG thumbnail w/ a border and return its aspect ratio.

    On error, returns None.
    """
    try:
        aspect = float(
            get_stdout(
                f'convert {file} -format "%[fx:w/h]" info:',
                shell=True,
            )
        )
        width = aspect * 500
        get_stdout(
            f"convert -resize {width}x500 -background white -alpha remove "
            f"-bordercolor black -border 15 {file} {file}.png",
            shell=True,
        )
    except:
        return None
    else:
        return aspect


if __name__ == "__main__":

    # Snakemake config (available automagically)
    config = snakemake.config  # type:ignore

    # User settings
    node_attr = config["dag"]["node_attr"]
    graph_attr = config["dag"]["graph_attr"]
    engine = config["dag"]["engine"]
    group_by_type = config["dag"]["group_by_type"]

    # Internal settings
    if group_by_type:
        subgraph_prefix = "cluster_"
        alpha = "33"
    else:
        subgraph_prefix = ""
        alpha = "ff"
    colors = dict(
        zenodo=f"#1f77b4{alpha}",  # blue
        other=f"#000000{alpha}",  # black
        script=f"#2ca02c{alpha}",  # green
        tex=f"#d62728{alpha}",  # red
        dataset=f"#9467bd{alpha}",  # purple
        figures="black",
        article=f"#ffffff{alpha}",  # white
        edge="black",
    )

    # Instantiate the graph
    dot = graphviz.Digraph(
        "dag",
        node_attr=node_attr,
        graph_attr=graph_attr,
        engine=engine,
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

    # Get permalink for file
    def file2url(file):
        return f"{config['git_url']}/blob/{config['git_sha']}/{file}"

    # Zenodo nodes
    with dot.subgraph(name=f"{subgraph_prefix}zenodo") as c:
        c.attr(
            color="black",
            fillcolor=colors["zenodo"],
            style="filled",
            label="https://doi.org/",
        )
        # Static
        for doi in get_dataset_dois(datasets, config["datasets"]):
            c.node(
                doi,
                URL=f"https://doi.org/{doi}",
                color="black" if group_by_type else colors["zenodo"],
                style="filled",
                fillcolor="white",
            )
        # Cache
        branch = config["git_branch"]
        cache = (
            config["cache"][branch]["zenodo"]
            or config["cache"][branch]["sandbox"]
        )
        if cache:
            c.node(
                cache,
                URL=f"https://doi.org/{cache}",
                color="black" if group_by_type else colors["zenodo"],
                style="filled",
                fillcolor="white",
            )

    # Dataset nodes
    with dot.subgraph(name=f"{subgraph_prefix}data") as c:
        c.attr(
            color="black",
            fillcolor=colors["dataset"],
            style="filled",
            label="src/data/",
        )
        for file in datasets:
            label = file.removeprefix("src/data/")
            c.node(
                file,
                label=label,
                color="black" if group_by_type else colors["dataset"],
                style="filled",
                fillcolor="white",
            )
            for url in get_dataset_dois([file], config["datasets"]):
                c.edge(doi, file)

    # Script nodes
    with dot.subgraph(name=f"{subgraph_prefix}scripts") as c:
        c.attr(
            color="black",
            fillcolor=colors["script"],
            style="filled",
            label="src/scripts/",
        )
        for file in scripts:
            label = file.removeprefix("src/scripts/")
            c.node(
                file,
                label=label,
                URL=file2url(file),
                color="black" if group_by_type else colors["script"],
                style="filled",
                fillcolor="white",
            )

    # Tex nodes
    with dot.subgraph(name=f"{subgraph_prefix}tex") as c:
        c.attr(
            color="black",
            fillcolor=colors["tex"],
            style="filled",
            label="src/tex/",
        )
        for file in texfiles:
            label = file.removeprefix("src/tex/")
            c.node(
                file,
                label=label,
                URL=file2url(file),
                color="black" if group_by_type else colors["tex"],
                style="filled",
                fillcolor="white",
            )

    # Figure nodes
    with dot.subgraph(name=f"{subgraph_prefix}figures") as c:
        c.attr(
            color="black",
            fillcolor="white",
            style="filled",
            label="src/tex/figures/",
        )
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
                    color="black" if group_by_type else colors["figures"],
                    style="filled",
                    fillcolor="white",
                )
            else:
                # Couldn't convert to PNG; display text instead
                label = file.removeprefix("src/data/")
                c.node(
                    file,
                    label=label,
                    color="black" if group_by_type else colors["figures"],
                    style="filled",
                    fillcolor="white",
                )

    # Other nodes
    with dot.subgraph(name=f"{subgraph_prefix}others") as c:
        c.attr(
            color="black", fillcolor=colors["other"], style="filled", label="/"
        )
        for file in others:
            c.node(
                file,
                URL=file2url(file),
                color="black" if group_by_type else colors["other"],
                style="filled",
                fillcolor="white",
            )

    # Article node
    with dot.subgraph(name=f"{subgraph_prefix}article") as c:
        c.attr(
            color="black",
            fillcolor=colors["article"],
            style="filled",
            label="ms.pdf",
        )
        c.node(
            config["ms_pdf"],
            label="",
            image=str(
                paths.showyourwork().resources / "img" / "article-thumb.png"
            ),
            penwidth="0",
            fixedsize="true",
            imagescale="false",
            width="1.15",
            height="1.5",
            color="black" if group_by_type else colors["article"],
            style="filled",
            fillcolor="white",
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
