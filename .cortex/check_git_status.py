from pathlib import Path
import subprocess
import glob


# Paths
ROOT = Path(__file__).parents[1].absolute()


def check_git_status(FIGURE_SCRIPTS=None):

    # By default, check all figure scripts
    if FIGURE_SCRIPTS is None:
        FIGURE_SCRIPTS = glob.glob(str(Path("figures") / "*.py"))

    # Check for figure scripts with uncommited changes
    for file in FIGURE_SCRIPTS:

        ctx_file = (
            ROOT / ".cortex" / "data" / "{}.cortex".format(Path(file).name)
        )

        # Check if the file is version controlled
        try:
            status = (
                subprocess.check_output(
                    ["git", "ls-files", "--error-unmatch", file], cwd=ROOT
                )
                .decode()
                .replace("\n", "")
            )
        except Exception as e:
            # File is not version controlled
            with open(ctx_file, "w") as f:
                print("ctxFigureScriptNotVersionControlled", file=f)
        else:
            # File is version controlled

            # Check if the file has uncommitted changes
            try:
                status = (
                    subprocess.check_output(
                        ["git", "status", "-s", file], cwd=ROOT
                    )
                    .decode()
                    .replace("\n", "")
                )
                if len(status):
                    raise Exception("Uncommitted changes!")
            except Exception as e:
                # Assume uncommited changes
                with open(ctx_file, "w") as f:
                    print("ctxFigureScriptHasUncommittedChanges", file=f)
            else:
                # File is good!
                with open(ctx_file, "w") as f:
                    print("ctxFigureScriptUpToDate", file=f)
