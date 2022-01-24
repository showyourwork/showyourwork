from utils import paths
import snakemake
import hashlib
import time


# Random unique flag file name
# Gets deleted on startup
blocking_flag_file = (
    (paths.temptemp / hashlib.sha256(str(time.time()).encode()).hexdigest())
    .relative_to(paths.user)
    .as_posix()
)


checkpoint blocking:
    """
    Dummy checkpoint to allow us to query the DAG before running other rules.

    """
    message:
        "Building DAG..."
    priority:
        snakemake.jobs.Job.HIGHEST_PRIORITY
    output:
        touch(blocking_flag_file)