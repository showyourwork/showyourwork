from showyourwork import paths
from showyourwork.config import get_upstream_dependencies
from showyourwork.patches import get_snakemake_variable
from showyourwork.zenodo import Zenodo, get_dataset_urls
import snakemake
from collections import defaultdict


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