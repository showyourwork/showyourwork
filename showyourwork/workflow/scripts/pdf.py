"""
Compiles the article manuscript into a PDF.

"""
from showyourwork import paths
from showyourwork.tex import compile_tex
from showyourwork.zenodo import get_dataset_urls
import sys
import shutil
from pathlib import Path
from jinja2 import Environment, BaseLoader


if __name__ == "__main__":

    # Snakemake config (available automagically)
    config = snakemake.config  # type:ignore

    # Get the Zenodo/Zenodo Sandbox cache public url
    branch = config["git_branch"]
    cache_zenodo_doi = config["cache"][branch]["zenodo"]
    cache_sandbox_doi = config["cache"][branch]["sandbox"]
    if cache_zenodo_doi or cache_sandbox_doi:
        if cache_zenodo_doi:
            record_id = cache_zenodo_doi.split("zenodo.")[1]
            zenodo_cache_url = f"https://zenodo.org/record/{record_id}"
        else:
            record_id = cache_sandbox_doi.split("zenodo.")[1]
            zenodo_cache_url = f"https://sandbox.zenodo.org/record/{record_id}"
    else:
        zenodo_cache_url = None

    # Do a final pass through the DAG to identify additional
    # dataset dependencies of each figure, and get their URLs
    # so we can add all the necessary margin icons
    recursive_dependencies = config["dag_dependencies_recursive"]
    for label, value in config["tree"]["figures"].items():

        # Recursively check figure dependencies for additional datasets
        datasets = value["datasets"]
        dependencies = value["dependencies"]
        upstream = set()
        for dep in dependencies:
            if dep in recursive_dependencies.keys():
                upstream |= set(recursive_dependencies[dep])
        upstream = list(upstream)

        # Add those datasets back to the config
        datasets = list(
            set(datasets) | set(get_dataset_urls(upstream, config["datasets"]))
        )
        config["tree"]["figures"][label]["datasets"] = datasets

        # If any of the upstream dependencies is a cached file,
        # we'll add a cache icon to the figure
        if set(upstream) & set(config["cached_deps"]):
            config["tree"]["figures"][label]["cached"] = True
        else:
            config["tree"]["figures"][label]["cached"] = False

    # Gather the figure script/dataset info so we can access it on the TeX side
    config["labels"] = {}
    for label, value in config["tree"]["figures"].items():
        # Figure script URL
        script = value["script"]
        if script is not None:
            config["labels"][f"{label}_script"] = script

        # Dataset URLs. Note that a max of 3 datasets will be displayed
        datasets = value["datasets"]
        for dataset, number in zip(datasets, ["One", "Two", "Three"]):
            config["labels"][f"{label}_dataset{number}"] = dataset

        # Cached output URL
        if value["cached"]:
            config["labels"][f"{label}_cache"] = zenodo_cache_url

    # Metadata file jinja template
    TEMPLATE = r"""
    ((* if github_actions *))
    \OnGithubActionstrue
    ((* else *))
    \OnGithubActionsfalse
    ((* endif *))
    \def\syw@url{((- git_url -))}
    \def\syw@sha{((- git_sha -))}
    \def\syw@runid{((- github_runid -))}

    ((* for key, value in labels.items() *))
    \addvalue{((- key -))}{((- value -))}
    ((* endfor *))

    % Check if the Git tag is set (i.e. is it an empty string)
    ((* if sha_tag_header != "" *))
    \newcommand{\gitHeader}{((- sha_tag_header -))}
    ((* endif *))
    """

    # Custom jinja environment for LaTeX
    ENV = Environment(
        block_start_string="((*",
        block_end_string="*))",
        variable_start_string="((-",
        variable_end_string="-))",
        comment_start_string="((=",
        comment_end_string="=))",
        trim_blocks=True,
        autoescape=False,
        loader=BaseLoader(),
    )

    # Generate the stylesheet metadata file
    with open(str(Path(snakemake.config["stylesheet_meta_file"])), "w") as f:
        meta = ENV.from_string(TEMPLATE).render(**snakemake.config)
        print(meta, file=f)

    # Build the paper
    compile_tex(
        snakemake.config,
        output_dir=paths.user().compile,
        stylesheet=paths.showyourwork().resources / "styles" / "build.tex",
    )

    # Copy the PDF to the user dir
    shutil.copy(
        str(paths.user().compile / (snakemake.config["ms_name"] + ".pdf")),
        str(Path(snakemake.config["ms_pdf"])),
    )
