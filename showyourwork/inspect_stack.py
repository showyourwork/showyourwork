"""
Implements functions that inspect the stack to infer things about the caller.

Specifically, implements `get_snakemake_variable` to access variables defined
within the Snakemake workflow so we can customize some showyourwork 
functionality.

This functionality is pretty hacky and is almost certainly not future-proof, 
so if we change the Snakemake version (currently pinned at 16.5.5), we may have 
to update the code in this file.

"""
from . import exceptions
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