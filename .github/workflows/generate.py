import subprocess
import jinja2
import shutil
import os
from pathlib import Path
from datetime import date

# Enviornment variables (set in .yml files)
TARGET_REPOSITORY = os.environ["SLUG"].split("/")[-1]
ACCESS_TOKEN = os.environ["ACCESS_TOKEN"]
REPO_ACTIVE = os.environ["REPO_ACTIVE"] == "true"
MANUSCRIPT_TYPE = os.environ["MANUSCRIPT_TYPE"]
RUN_TESTS = os.environ["RUN_TESTS"] == "true"
AUTHOR = os.environ["AUTHOR"]
SLUG = os.environ["SLUG"]
README_MESSAGE = os.environ["README_MESSAGE"]
SHA = os.environ["SHA"]

# Showyourwork version. If there are any commits on top of the latest tag,
# use the exact SHA for the current commit on rodluger/showyourwork.
# Otherwise, use the exact tag name
try:
    tag_and_commits = (
        subprocess.check_output(["git", "describe"], stderr=subprocess.DEVNULL)
        .decode()
        .replace("\n", "")
    )
except Exception as e:
    print(e)
    tag_and_commits = ""
try:
    tag = (
        subprocess.check_output(
            ["git", "describe", "--abbrev=0"], stderr=subprocess.DEVNULL
        )
        .decode()
        .replace("\n", "")
    )
except Exception as e:
    print(e)
    tag = ""
if tag_and_commits == tag:
    version = tag
else:
    version = SHA

# Jinja keywords
kwargs = dict(
    version=version,
    manuscript_type=MANUSCRIPT_TYPE,
    year=date.today().year,
    author=AUTHOR,
    readme_message=README_MESSAGE,
    slug=SLUG,
    repo_active=REPO_ACTIVE,
    run_tests=RUN_TESTS,
)

# Create the new repo
shutil.copytree(
    Path("showyourwork") / "templates" / "repository", TARGET_REPOSITORY
)

# Render the templates
env = jinja2.Environment(trim_blocks=True, loader=jinja2.FileSystemLoader("."))
for file in Path(TARGET_REPOSITORY).glob("**/*"):
    file = str(file)
    if os.path.isfile(file):
        sty = env.get_template(file).render(**kwargs)
        with open(file, "w") as f:
            print(sty, file=f)

# Create repo and force push to target
subprocess.check_call(["git", "init"], cwd=TARGET_REPOSITORY)
subprocess.check_call(
    ["git", "checkout", "--orphan", "main"], cwd=TARGET_REPOSITORY
)
subprocess.check_call(["git", "add", "."], cwd=TARGET_REPOSITORY)
subprocess.check_call(
    [
        "git",
        "-c",
        "user.name='showyourwork'",
        "-c",
        "user.email='showyourwork",
        "commit",
        "-m",
        "'auto commit from showyourwork'",
    ],
    cwd=TARGET_REPOSITORY,
)
subprocess.check_call(
    [
        "git",
        "push",
        "--force",
        f"https://x-access-token:{ACCESS_TOKEN}@github.com/{SLUG}",
        "main",
    ],
    cwd=TARGET_REPOSITORY,
    stdout=subprocess.DEVNULL,
    stderr=subprocess.STDOUT,
)
