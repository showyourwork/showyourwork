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

        # Temporary metadata file
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
            error_code = "ctxFigureScriptNotVersionControlled"
        else:
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
                error_code = "ctxFigureScriptHasUncommittedChanges"
            else:
                # File is good!
                error_code = "ctxFigureScriptUpToDate"

        if ctx_file.exists():
            with open(ctx_file, "r") as f:
                ctx_file_contents = f.read().replace("\n", "")
        else:
            ctx_file_contents = ""
        if ctx_file_contents != error_code:
            with open(ctx_file, "w") as f:
                print(error_code, file=f)
