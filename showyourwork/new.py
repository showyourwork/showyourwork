import subprocess
from pathlib import Path

HERE = Path(__file__).absolute().parents[0]


def new():
    subprocess.check_call(["cookiecutter", HERE / "cookiecutter-showyourwork"])
