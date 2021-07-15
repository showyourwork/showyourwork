import subprocess
from pathlib import Path
from glob import glob
import yaml
import os
import shutil
import argparse


ROOT = Path(
    os.getenv("GITHUB_WORKSPACE", Path(__file__).parents[2].absolute())
)

ANCIENT = "200001010000.00"  # Jan 1 2000 00:00:00


def get_figure_extensions():
    """
    Get user-defined figure extensions.

    """
    with open(ROOT / "config.yml", "r") as f:
        config = yaml.safe_load(f)
    return config.get(
        "figure_extensions", ["pdf", "png", "epsfont", "jpg", "jpeg", "gif"]
    )


def get_modified_files(commit="HEAD^"):
    """
    Return all files that changed since `commit`.

    """
    return [
        file
        for file in (
            subprocess.check_output(
                ["git", "diff", "HEAD", commit, "--name-only"]
            )
            .decode()
            .split("\n")
        )
        if len(file)
    ]


def restore_cache():
    """
    Restore cached figure files and test results.

    """
    # Copy the cached figure files
    for ext in get_figure_extensions():
        for file in glob(
            str(ROOT / ".cache" / "figures" / "*.{}".format(ext))
        ):
            shutil.copy(file, ROOT / "figures")
            print(
                "Restored file figures/{} from cache.".format(Path(file).name)
            )

    # Cache test results
    for file in glob(str(ROOT / ".cache" / "tests" / "*.status")):
        shutil.copy(file, ROOT / "tests")
        print("Restored file tests/{} from cache.".format(Path(file).name))

    # Get the commit when the files were cached
    try:
        with open(ROOT / ".cache" / "commit", "r") as f:
            commit = f.readlines()[0].replace("\n", "")
    except FileNotFoundError:
        return

    # Get all files modified since that commit
    files = get_modified_files(commit)

    # If a script has not changed since the commit at which its
    # output was cached, mark it as ancient to prevent Snakemake
    # from re-running it.
    for script in glob(str(ROOT / "figures" / "*.py")):
        script_rel = str(Path("figures") / script.name)
        if script_rel not in files:
            subprocess.check_output(
                ["touch", "-a", "-m", "-t", ANCIENT, script]
            )
            print("File `{}` marked as ANCIENT.".format(script_rel))

    for script in glob(str(ROOT / "tests" / "test_*.py")):
        script_rel = str(Path("tests") / script.name)
        if script_rel not in files:
            subprocess.check_output(
                ["touch", "-a", "-m", "-t", ANCIENT, script]
            )
            print("File `{}` marked as ANCIENT.".format(script_rel))


def update_cache():
    """
    Cache all figure files and test results.

    """
    # Create the cache folder
    if not Path(ROOT / ".cache").exists():
        os.mkdir(ROOT / ".cache")

    # Cache figure files
    for ext in get_figure_extensions():
        for file in glob(str(ROOT / "figures" / "*.{}".format(ext))):
            shutil.copy(file, ROOT / ".cache" / "figures")
            print("Cached file figures/{}.".format(Path(file).name))

    # Cache test results
    for file in glob(str(ROOT / "tests" / "*.status")):
        shutil.copy(file, ROOT / ".cache" / "tests")
        print("Cached file tests/{}.".format(Path(file).name))

    # Store the current commit
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode()
    with open(ROOT / ".cache" / "commit", "w") as f:
        print(commit, file=f)


if __name__ == "__main__":

    # Command line options
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--restore", action="store_true")
    parser.add_argument("-u", "--update", action="store_true")
    args = parser.parse_args()
    if args.restore:
        restore_cache()
    elif args.update:
        update_cache()
    else:
        raise ValueError("Please provide either --restore or --update.")
