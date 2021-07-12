import subprocess
from pathlib import Path
import re
import glob
import os
import json


# Paths
ROOT = Path(__file__).parents[1].absolute()


def infer_figure_scripts():

    with open(ROOT / "tex" / "dryrun.tex", "w") as f:
        print("\def\ctxDryRun{}\input{ms}", file=f)
    subprocess.check_output(
        ["tectonic", "--keep-logs", "dryrun.tex"], cwd=str(ROOT / "tex")
    )
    with open(ROOT / "tex" / "cortex-dryrun.log", "r") as f:
        log = f.readlines()

    scripts = {}
    for line in log:
        match = re.match("FIGURE: `(.*?)` --> `(.*?)`.", line)
        if match:
            script, file = match.groups()
            script = str(Path("figures") / script)
            file = str(Path("figures") / file)
            if script in scripts:
                scripts[script] += [file]
            else:
                scripts[script] = [file]
    for file in glob.glob(str(ROOT / "tex" / "dryrun.*")):
        os.remove(file)

    with open(ROOT / ".cortex" / "data" / "figure_scripts.json", "w") as f:
        json.dump(scripts, f)
