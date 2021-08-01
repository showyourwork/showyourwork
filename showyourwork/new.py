import subprocess


def new(version="latest"):
    if version == "latest":
        subprocess.check_call(
            ["cookiecutter", "gh:rodluger/cookiecutter-showyourwork"]
        )
    else:
        if not version.startswith("v"):
            version = "v" + version
        subprocess.check_call(
            [
                "cookiecutter",
                f"https://github.com/rodluger/cookiecutter-showyourwork/archive/refs/tags/{version}.zip",
            ]
        )
