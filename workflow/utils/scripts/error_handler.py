import platform
import os
import shutil
import sys
from pathlib import Path


WORKFLOW = Path(__file__).absolute().parents[2]
sys.path.insert(1, str(WORKFLOW))
from utils import paths


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
    for file in [paths.logs / "preprocess.log", paths.logs / "compile.log"]:
        if file.exists():
            with open(file, "r") as f:
                line = f.readline().replace("\n", "")
            if line.endswith(".snakemake.log") and Path(line).exists():
                shutil.copy(Path(line), file)


def print_errors(errorcode):

    # Get showyourwork errors
    syw_error_file = paths.logs / "showyourwork_errors.log"
    if syw_error_file.exists():
        with open(syw_error_file, "r") as f:
            syw_errors = f.read()
            if syw_errors.endswith("\n"):
                syw_errors = syw_errors[:-1]
    else:
        syw_errors = ""

    # Get snakemake errors
    sm_error_file = paths.logs / "snakemake_errors.log"
    if sm_error_file.exists():
        with open(sm_error_file, "r") as f:
            sm_errors = f.read()
            if sm_errors.endswith("\n"):
                sm_errors = sm_errors[:-1]
    else:
        sm_errors = ""

    # Print the errors. Note that showyourwork errors
    # take precedence.
    if len(syw_errors):
        # We already printed the error when it was raised
        sys.exit(errorcode)
    elif len(sm_errors):
        print(colorize(sm_errors))
        sys.exit(errorcode)
    else:
        sys.exit(int(errorcode))


if __name__ == "__main__":
    assert len(sys.argv) == 2, "Incorrect number of args to `print_errors.py`."
    errorcode = int(sys.argv[1])
    copy_logfiles()
    print_errors(errorcode)