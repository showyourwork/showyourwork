import graphviz
from PyPDF2 import PdfFileReader
from pathlib import Path
import json


# Node border colors
colors = dict(
    zenodo="#1f77b4",
    analysis="black",
    script="black",
    figure="black",
    dataset="#1f77b4",
    other="black",
    edge="black",
)

# Get figure metadata
dependencies = snakemake.config.get("dependencies", {})
zenodo = snakemake.config.get("zenodo", {})
sandbox = snakemake.config.get("zenodo_sandbox", {})
FIGURES = snakemake.params.get("FIGURES", "src/figures")
with open(".showyourwork/scripts.json", "r") as f:
    script_info = json.load(f)

# Normalize the entries the user provides them as an OrderedDict
# in the YAML file (i.e., as list `-` entries)
if type(zenodo) is list:
    zz = {}
    for z in zenodo:
        zz.update(dict(z))
    zenodo = zz
if type(sandbox) is list:
    zz = {}
    for z in sandbox:
        zz.update(dict(z))
    sandbox = zz

# Instantiate the graph
dot = graphviz.Digraph("dag", node_attr={"shape": "box", "penwidth": "2", "width": "1"})

# Figures
png_files = []
for figure in script_info.get("figures", []):

    # Get metadata for this figure
    script = script_info["figures"][figure]["script"]
    files = script_info["figures"][figure]["files"]

    # Figure node
    dot.node(script, label=script.replace("src/", ""), color=colors["script"])

    # Loop through figure dependencies
    for dependency in dependencies.get(script, []):
        if not dependency.endswith(".zenodo"):

            # Find the relevant Zenodo or Zenodo tarball entry
            # (if this is a Zenodo-hosted dependency)
            zenodo_entry = zenodo.get(dependency, {})
            if not zenodo_entry:
                for entry in zenodo:
                    if dependency in zenodo[entry].get("contents", []):
                        zenodo_entry = zenodo[entry]
                        break
            sandbox_entry = sandbox.get(dependency, {})
            if not sandbox_entry:
                for entry in sandbox:
                    if dependency in sandbox[entry].get("contents", []):
                        sandbox_entry = sandbox[entry]
                        break

            # Dependency node
            if zenodo_entry or sandbox_entry:
                dot.node(
                    dependency,
                    label=dependency.replace("src/", ""),
                    color=colors["dataset"],
                    shape="box3d",
                )
            else:
                dot.node(
                    dependency,
                    label=dependency.replace("src/", ""),
                    color=colors["other"],
                    style="rounded",
                )

            # Zenodo metadata for this dependency
            for zenodo_info, zenodo_stem in zip(
                [zenodo_entry, sandbox_entry],
                ["10.5281", "10.5072"],
            ):
                if zenodo_info:

                    # Zenodo node
                    zenodo_id = zenodo_info["id"]
                    doi = f"{zenodo_stem}/zenodo.{zenodo_id}"
                    dot.node(
                        doi,
                        label=doi.replace("src/", ""),
                        color=colors["zenodo"],
                        shape="cylinder",
                        height="0.75",
                    )
                    dot.edge(doi, dependency, color=colors["edge"])

                    # Parent script node
                    zenodo_script = zenodo_info.get("script")
                    if zenodo_script:
                        dot.node(
                            zenodo_script,
                            label=zenodo_script.replace("src/", ""),
                            color=colors["analysis"],
                            style="rounded",
                        )
                        dot.edge(
                            zenodo_script, doi, color=colors["edge"], style="dashed"
                        )

                        breakpoint()

                        dot.edge(
                            zenodo_script,
                            dependency,
                            color=colors["edge"],
                            style="dashed",
                        )

            # Connect stuff
            for file in files:
                dot.edge(dependency, file, color=colors["edge"])

    # Figure images
    for file in files:

        # DOT can't render PDF images, so we'll convert them to PNG later
        file_png = file[: -len(Path(file).suffix)] + ".png"

        # If the file is a PDF, let's infer its aspect ratio
        # TODO: We could do this for other figure types as well... PR anyone?
        if Path(file).suffix.lower() == ".pdf":
            with open(file, "rb") as f:
                w, h = PdfFileReader(f).getPage(0).mediaBox[2:]
            aspect = float(w / h)
        else:
            aspect = 2

        dot.node(
            file,
            label="",
            image=file_png,
            penwidth="4",
            fixedsize="true",
            imagescale="false",
            width=f"{aspect}",
            height="1",
            color=colors["figure"],
        )
        dot.edge(script, file, color=colors["edge"])

        # Invisible edge (to force same rank for all figures)
        dot.edge(file, "anchor", color="#00000000")

# Invisble anchor
dot.node("anchor", label="", height="0", color="#00000000")

# Save the graph
dot.save(directory=FIGURES)