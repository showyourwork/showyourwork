from showyourwork import paths
from showyourwork.patches import get_snakemake_variable
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

        # Assemble input-output information; add it to the global config
        rule_dependencies = defaultdict(set)
        for rule in dag.rules:
            outputs = [str(file) for file in rule.output]
            inputs = [
                str(file) for file in rule.input 
                if type(file) is snakemake.io._IOFile
            ]
            for file in outputs:
                rule_dependencies[file] |= set(inputs)
        for key in rule_dependencies.keys():
            rule_dependencies[key] = list(rule_dependencies[key])
        snakemake.workflow.config["rule_dependencies"] = dict(rule_dependencies)

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