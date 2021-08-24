from .constants import *
from pathlib import Path
from glob import glob as _glob
import shutil
import subprocess
import json
import jinja2
import sys
import re


__all__ = ["glob", "save_json", "make_pdf", "check_repo"]


def glob(pathname, **kwargs):
    """
    Same as `glob.glob`, but works on strings and `pathlib.Path` instances.

    """
    return _glob(str(pathname), **kwargs)


def make_pdf(
    tmpdir=TEMP / "tex",
    publish=True,
    verbose=False,
    tectonic_args=["--keep-logs", "--keep-intermediates"],
    **jinja_kwargs,
):
    """
    Compile the PDF using tectonic.

    """
    # Custom Jinja environment for LaTeX
    env = jinja2.Environment(
        block_start_string="((*",
        block_end_string="*))",
        variable_start_string="((-",
        variable_end_string="-))",
        comment_start_string="((=",
        comment_end_string="=))",
        trim_blocks=True,
        autoescape=False,
        loader=jinja2.FileSystemLoader(ROOT / "templates"),
    )

    # The style sheet template
    style = "showyourwork.sty"

    # Generate the stylesheet
    with open(USER / "tex" / f"{style}", "w") as f:
        sty = env.get_template(style).render(**jinja_kwargs)
        print(sty, file=f)

    # Copy LaTeX stylesheets & helper files to tmpdir
    for file in glob(ROOT / "tex" / "*"):
        shutil.copy(file, tmpdir)

    # Copy user files to tmpdir
    for file in glob(USER / "tex" / "*"):
        shutil.copy(file, tmpdir)

    # Include the showyourwork package right after `\documentclass` call
    with open(tmpdir / "ms.tex", "r") as f:
        lines = f.readlines()
    for idx, line in enumerate(lines):
        if line.startswith(r"\documentclass"):
            lines = (
                lines[: idx + 1]
                + [r"\usepackage{showyourwork}" + "\n"]
                + lines[idx + 1 :]
            )
            break
    else:
        raise ValueError(r"Missing `\documentclass` in file `tex/ms.tex`.")
    with open(tmpdir / "ms.tex", "w") as f:
        f.writelines(lines)

    # Verbosity
    if verbose:
        tectonic_args += ["--print"]
    else:
        tectonic_args += ["--chatter", "minimal"]

    # Build
    subprocess.check_call(
        ["tectonic"] + tectonic_args + ["ms.tex"], cwd=tmpdir
    )

    # Copy PDF to user directory
    if publish:
        shutil.copy2(tmpdir / "ms.pdf", USER)


def save_json(contents, file):
    """
    Saves a dictionary in JSON format to disk only if the contents
    of the file would be changed (or if the file doesn't exist).

    """
    file = Path(file)
    current = json.dumps(contents, indent=4)
    if file.exists():
        try:
            with open(file, "r") as f:
                existing = json.load(f)
        except json.JSONDecodeError as e:
            existing = None
        if existing != current:
            with open(file, "w") as f:
                print(current, file=f)
    else:
        file.parents[0].mkdir(exist_ok=True, parents=True)
        with open(file, "w") as f:
            print(current, file=f)


def check_repo():
    if not (Path("tex") / "ms.tex").exists():
        print("ERROR: Missing manuscript file `tex/ms.tex`.")
        sys.exit(1)
    if not Path("figures").exists():
        print("ERROR: Missing figure scripts directory `figures`.")
        sys.exit(1)
    if not Path("environment.yml").exists():
        print("ERROR: Missing conda environment file `environment.yml`.")
        sys.exit(1)
