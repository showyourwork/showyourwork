import shutil
from pathlib import Path

from . import exceptions, paths
from .logging import get_logger
from .subproc import get_stdout

__all__ = ["compile_tex"]


def compile_tex(config, output_dir=None, args=[], stylesheet=None):
    """
    Compile the TeX document using ``tectonic``.

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
        f'"{output_dir}"',
    ]

    def callback(code, stdout, stderr):
        # Copy over the tectonic log file
        logfile = Path(output_dir) / "{}.log".format(config["ms_name"])
        if logfile.exists():
            shutil.copy(logfile, paths.user().logs / "tectonic.log")
            logfile = paths.user().logs / "tectonic.log"
        else:
            logfile = None

        if code > 0:
            # Log the error
            logger = get_logger()
            logger.error(stderr)
            # Raise the exception
            raise exceptions.TectonicError(logfile)

    # Activate the showyourwork conda environment
    conda_activate = open(paths.user().flags / "SYW__CONDA", "r").read()

    get_stdout(
        f"{conda_activate} tectonic {' '.join(args)} {' '.join(force_args)} {' '.join(config['user_args'])} \"{paths.user().repo / config['ms_tex']}\"",
        shell=True,
        callback=callback,
    )
