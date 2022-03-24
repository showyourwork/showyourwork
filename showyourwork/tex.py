from . import paths, exceptions
from .logging import get_logger
from pathlib import Path
import subprocess
import shutil


__all__ = ["compile_tex"]


def compile_tex(config, output_dir=None, args=[], stylesheet=None):
    """
    Compile the TeX document using `tectonic`.

    """
    if output_dir is None:
        output_dir = paths.user().compile

    # Copy over TeX auxiliaries
    for file in config["tex_files_in"]:
        src = paths.user().repo / file
        dst = paths.user().tex / src.name
        if not dst.exists():
            shutil.copy(str(src), str(dst))

    # Copy over the actual stylesheet
    if stylesheet is not None:
        shutil.copy(
            str(stylesheet), str(paths.user().tex / "showyourwork.tex")
        )

    # Run tectonic
    force_args = [
        "--chatter",
        "minimal",
        "--keep-logs",
        "--keep-intermediates",
        "-o",
        str(output_dir),
    ]
    result = subprocess.run(
        ["tectonic"]
        + args
        + force_args
        + [paths.user().repo / config["ms_tex"]],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    # Copy over the tectonic log file
    logfile = Path(output_dir) / "{}.log".format(config["ms_name"])
    if logfile.exists():
        shutil.copy(logfile, paths.user().logs / "tectonic.log")
        logfile = paths.user().logs / "tectonic.log"
    else:
        logfile = None

    # Process errors
    if result.returncode != 0:

        # Log the error
        logger = get_logger()
        logger.error(result.stderr.decode("utf-8"))

        # Raise the exception
        with exceptions.no_traceback():
            raise exceptions.TectonicError(logfile)