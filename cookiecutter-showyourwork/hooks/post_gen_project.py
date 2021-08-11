import subprocess

if "{{ cookiecutter.access_token }}" != "":

    # Create repo and force push to target
    subprocess.check_call(["git", "init"])
    subprocess.check_call(["git", "checkout", "--orphan", "main"])
    subprocess.check_call(["git", "add", "."])
    if "{{ cookiecutter.skip_ci }}" == "n":
        message = "auto commit from rodluger/showyourwork"
    else:
        message = "[skip ci] auto commit from rodluger/showyourwork"
    subprocess.check_call(
        [
            "git",
            "-c",
            "user.name=rodluger/showyourwork",
            "-c",
            "user.email=rodluger/showyourwork",
            "commit",
            "-m",
            message,
        ],
    )
    subprocess.check_call(
        [
            "git",
            "push",
            "--force",
            "https://x-access-token:{{ cookiecutter.access_token }}@github.com/{{ cookiecutter.slug }}",
            "main",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )
