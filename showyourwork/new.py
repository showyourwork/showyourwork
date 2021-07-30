import subprocess


def new():
    subprocess.check_call(
        ["cookiecutter", "gh:rodluger/cookiecutter-showyourwork"]
    )
