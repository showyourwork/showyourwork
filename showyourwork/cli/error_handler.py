from .. import paths
import platform
import os
import shutil
import sys
from pathlib import Path


def colorize(message, color="red"):
    colors = {
        "black": 30,
        "red": 31,
        "green": 32,
        "yellow": 33,
        "blue": 34,
        "magenta": 35,
        "cyan": 36,
        "white": 37,
    }
    if platform.system() == "Windows" or os.environ.get("TERM", "") == "dumb":
        return message
    else:
        return f"\033[{colors[color]}m{message}\033[0m"


def copy_logfiles():
    for file in [
        paths.user().logs / "preprocess.log",
        paths.user().logs / "compile.log",
    ]:
        if file.exists():
            with open(file, "r") as f:
                line = f.readline().replace("\n", "")
            if line.endswith(".snakemake.log") and Path(line).exists():
                shutil.copy(Path(line), file)


def print_errors(errorcode):

    # Get showyourwork errors
    syw_error_file = paths.user().logs / "showyourwork_errors.log"
    if syw_error_file.exists():
        with open(syw_error_file, "r") as f:
            syw_errors = f.read()
            if syw_errors.endswith("\n"):
                syw_errors = syw_errors[:-1]
    else:
        syw_errors = ""

    # Get snakemake errors
    sm_error_file = paths.user().logs / "snakemake_errors.log"
    if sm_error_file.exists():
        with open(sm_error_file, "r") as f:
            sm_errors = f.read()
            if sm_errors.endswith("\n"):
                sm_errors = sm_errors[:-1]
    else:
        sm_errors = ""

    # Print any errors and exit on failure.
    # Note that showyourwork errors take precedence.
    if len(syw_errors):
        # We already printed the error when it was raised
        sys.exit(errorcode)
    elif len(sm_errors):
        print(colorize(sm_errors))
        sys.exit(errorcode)
    elif int(errorcode) > 0:
        sys.exit(int(errorcode))


def error_handler(errorcode):
    copy_logfiles()
    print_errors(errorcode)