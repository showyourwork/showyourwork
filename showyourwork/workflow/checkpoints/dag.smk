from showyourwork import paths
import snakemake


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