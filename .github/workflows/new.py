from pathlib import Path
from cookiecutter.main import cookiecutter
from datetime import date
import subprocess
from packaging.version import parse as parse_version
import os
import sys


# Root of the `showyourwork` repository
ROOT = Path(__file__).absolute().parents[2]


def new(slug, repo_active="y", access_token="", skip_ci="n", run_tests="n"):

    # Parse user args
    assert (
        len(slug.split("/")) == 2
    ), "Argument `slug` must have the form `user/repo`."
    user, repo = slug.split("/")

    # Infer showyourwork version
    try:
        # Use the package version, if available and if release
        from showyourwork import __version__

        assert not parse_version(__version__).is_prerelease
        version = __version__
    except:
        try:
            # Use the latest commit SHA, if available
            version = (
                subprocess.check_output(
                    ["git", "rev-parse", "HEAD"],
                    stderr=subprocess.DEVNULL,
                    cwd=ROOT,
                )
                .decode()
                .replace("\n", "")
            )
        except Exception as e:
            # Fallback
            version = "main"

    # Coookiecutter config
    extra_context = {
        "slug": slug,
        "user": user,
        "repo": repo,
        "year": date.today().year,
        "version": version,
        "repo_active": repo_active,
        "access_token": access_token,
        "skip_ci": skip_ci,
        "run_tests": run_tests,
    }

    # Render
    cookiecutter(
        str(ROOT / "cookiecutter-showyourwork"),
        extra_context=extra_context,
        no_input=True,
    )


if __name__ == "__main__":

    if os.environ.get("CI", "false") == "true":

        # Running on GitHub Actions
        slug = os.environ["SLUG"]
        repo_active = os.environ["REPO_ACTIVE"]
        access_token = os.environ["ACCESS_TOKEN"]
        skip_ci = os.environ["SKIP_CI"]
        run_tests = os.environ["RUN_TESTS"]
        new(
            slug,
            repo_active=repo_active,
            access_token=access_token,
            skip_ci=skip_ci,
            run_tests=run_tests,
        )

    else:

        # Running locally from the command line
        assert (
            len(sys.argv) == 2
        ), "Please provide a single argument corresponding to the repo slug."
        new(sys.argv[1])
