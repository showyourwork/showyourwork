from .showyourwork_version import __version__
import subprocess
import jinja2
import shutil
import os
from pathlib import Path
from datetime import date
from packaging.version import parse as parse_version


ROOT = Path(__file__).absolute().parents[1]


def new(slug):
    # Arg check
    assert (
        len(slug.split("/")) == 2
    ), "Argument `slug` must have the form `user/repo`."

    # Infer showyourwork version
    if parse_version(__version__).is_prerelease:
        try:
            tag_and_commits = (
                subprocess.check_output(
                    ["git", "describe"], stderr=subprocess.DEVNULL, cwd=ROOT
                )
                .decode()
                .replace("\n", "")
            )
        except Exception as e:
            tag_and_commits = "tag_and_commits"
        try:
            tag = (
                subprocess.check_output(
                    ["git", "describe", "--abbrev=0"],
                    stderr=subprocess.DEVNULL,
                    cwd=ROOT,
                )
                .decode()
                .replace("\n", "")
            )
        except Exception as e:
            tag = "tag"
        if tag_and_commits == tag:
            version = tag
        else:
            version = "main"  # fallback!
    else:
        version = __version__

    # Jinja keywords
    author = slug.split("/")[0]
    kwargs = dict(
        version=version,
        manuscript_type="minimal",
        year=date.today().year,
        author=author,
        readme_message="",
        slug=slug,
        repo_active=True,
        run_tests=False,
    )

    # Create the new repo
    repo = slug.split("/")[-1]
    shutil.copytree(ROOT / "showyourwork-template", repo)

    # Render the templates
    env = jinja2.Environment(
        trim_blocks=True, loader=jinja2.FileSystemLoader(".")
    )
    for file in Path(repo).glob("**/*"):
        file = str(file)
        if os.path.isfile(file):
            sty = env.get_template(file).render(**kwargs)
            with open(file, "w") as f:
                print(sty, file=f)

    #
    print(f"Created article directory `{repo}`.")
