import inspect
import sys

try:
    import snakemake
except ModuleNotFoundError:
    snakemake = None


def get_snakemake_variable(name, default=None):
    """
    Infer the value of a variable within snakemake.

    This is extremely hacky.
    """
    levels = inspect.stack()
    for level in levels:
        value = level.frame.f_locals.get(name, None)
        if value is not None:
            return value
    return default


def is_make_clean():
    """
    Returns True if the current workflow was triggered by `make clean`

    """
    if "--delete-all-output" in sys.argv:
        return True
    else:
        return False


def is_make_main():
    """
    Returns True if the current workflow is the main PDF build.

    """
    dag = get_snakemake_variable("dag")
    targetjobs = [j.name for j in dag.targetjobs]
    if "syw__main" in targetjobs or "syw__compile" in targetjobs:
        return True
    else:
        return False