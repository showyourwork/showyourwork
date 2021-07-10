from pathlib import Path
import subprocess
import glob


# Paths
ROOT = Path(__file__).parents[1].absolute()

# Check for figure scripts with uncommited changes
for file in glob.glob(str(ROOT / "figures" / "*.py")):

    # Check if the file is version controlled
    try:
        status = (
            subprocess.check_output(
                ["git", "ls-files", "--error-unmatch", file]
            )
            .decode()
            .replace("\n", "")
        )
    except Exception as e:
        # File is not version controlled
        with open("{}.cortex".format(file), "w") as f:
            print("ctxFigureScriptNotVersionControlled", file=f)
    else:
        # File is version controlled

        # Check if the file has uncommitted changes
        try:
            status = (
                subprocess.check_output(["git", "status", "-s", file])
                .decode()
                .replace("\n", "")
            )
            if len(status):
                raise Exception("Uncommitted changes!")
        except Exception as e:
            # Assume uncommited changes
            with open("{}.cortex".format(file), "w") as f:
                print("ctxFigureScriptHasUncommittedChanges", file=f)
        else:
            # File is good!
            with open("{}.cortex".format(file), "w") as f:
                print("ctxFigureScriptUpToDate", file=f)
