import os
import shutil

from ... import paths
from .run_snakemake import run_snakemake


# To handle WindowsError: [Error 5] Access is denied:  https://stackoverflow.com/a/2656405
def onerror(func, path, exc_info):
    """
    Error handler for ``shutil.rmtree``.

    If the error is due to an access error (read only file)
    it attempts to add write permission and then retries.

    If the error is for another reason it re-raises the error.

    Usage : ``shutil.rmtree(path, onerror=onerror)``
    """
    import stat

    # Is the error an access error?
    if not os.access(path, os.W_OK):
        os.chmod(path, stat.S_IWUSR)
        func(path)
    else:
        raise OSError(f"Could not delete {path}")


def clean(force, deep, snakemake_args=(), cores=1, conda_frontend="conda"):
    """Clean the article build.

    Args:
        force (bool): If True, forcefully delete files in output directories.
        deep (bool): If True, delete all temporary Snakemake and showyourwork
            directories.
        options (str, optional): Additional options to pass to Snakemake.

    """

    if (paths.user().repo / ".snakemake" / "incomplete").exists():
        shutil.rmtree(paths.user().repo / ".snakemake" / "incomplete", onerror=onerror)
    for file in ["build.smk", "prep.smk"]:
        snakefile = paths.showyourwork().workflow / file
        run_snakemake(
            snakefile.as_posix(),
            run_type="clean",
            cores=cores,
            conda_frontend=conda_frontend,
            extra_args=list(snakemake_args) + ["--delete-all-output"],
            check=False,
        )
    if (paths.user().repo / "arxiv.tar.gz").exists():
        (paths.user().repo / "arxiv.tar.gz").unlink()
    if paths.user().temp.exists():
        shutil.rmtree(paths.user().temp, onerror=onerror)
    if force:
        for root, dirs, files in os.walk(paths.user().data, topdown=False):
            for name in files:
                if name != ".gitignore":
                    os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        for root, dirs, files in os.walk(paths.user().figures, topdown=False):
            for name in files:
                if name != ".gitignore":
                    os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
    if deep:
        if (paths.user().repo / ".snakemake").exists():
            shutil.rmtree(paths.user().repo / ".snakemake", onerror=onerror)
        if (paths.user().repo / ".showyourwork").exists():
            shutil.rmtree(paths.user().repo / ".showyourwork", onerror=onerror)
