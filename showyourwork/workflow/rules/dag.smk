"""
Defines the checkpoint ``syw__dag`` to query the Snakemake DAG.

This checkpoint runs before all the other rules so we can programmatically
access the workflow graph and use information about it within the workflow
itself.

"""
from showyourwork import paths
from showyourwork.config import get_upstream_dependencies
from showyourwork.patches import get_snakemake_variable
from showyourwork.zenodo import Zenodo, get_dataset_urls
import snakemake
from collections import defaultdict

def infer_additional_figure_dependencies():
    """
    Infers additional figure dependencies directly from the Snakemake DAG.

    These are added to the config so we have a complete picture of all
    files upstream of the figures (for margin icon generation & DAG
    rendering).

    """
    # Snakemake config
    config = snakemake.workflow.config

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
        if set(dependencies) & set(config["cached_deps"]):
            config["tree"]["figures"][label]["cached"] = True
        elif set(upstream) & set(config["cached_deps"]):
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


def WORKFLOW_GRAPH(*args):
    """
    This is a dummy input function for the main rules of the build (``syw__pdf``
    and ``syw_arxiv``) which returns an empty list. As a side effect, it blocks
    until the checkpoint ``syw_dag`` is executed, at which point it queries the
    workflow dag and adds graph metadata to the global ``config``.

    """
    # Wait until the DAG has been built
    snakemake.workflow.checkpoints.syw__dag.get()

    # Get the workflow graph
    dag = get_snakemake_variable("dag", None)
    
    if dag is None:
        
        # Fail silently
        pass

    else:

        # Assemble input-output information
        dependencies = defaultdict(set)
        for job in dag.jobs:
            outputs = [str(file) for file in job.output]
            inputs = [
                str(file) for file in job.input 
                if type(file) is snakemake.io._IOFile
            ]
            for file in outputs:
                dependencies[file] |= set(inputs)
        for key in dependencies.keys():
            dependencies[key] = list(dependencies[key])
        dependencies = dict(dependencies)
        
        # Find recursive input-output dependencies
        recursive_dependencies = {}
        for output in dependencies.keys():
            recursive_dependencies[output] = get_upstream_dependencies(
                output, 
                dependencies
            )
        
        # Add to the global config
        config = snakemake.workflow.config
        config["dag_dependencies"] = dependencies
        config["dag_dependencies_recursive"] = recursive_dependencies

        # Infer figure deps from the DAG
        infer_additional_figure_dependencies()

    # Dummy output
    return []


checkpoint:
    """
    Dummy checkpoint to allow us to query the DAG before running other rules.

    """
    name:
        "syw__dag"
    message:
        "Building DAG..."
    priority:
        snakemake.jobs.Job.HIGHEST_PRIORITY
    output:
        touch(paths.user().flags / "SYW__DAG")