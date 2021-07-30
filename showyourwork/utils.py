from . import settings
from .constants import *
from pathlib import Path
from glob import glob as _glob
import shutil
import subprocess
import json
import jinja2


__all__ = ["glob", "save_json", "make_pdf"]


def glob(pathname, **kwargs):
    """
    Same as `glob.glob`, but works on strings and `pathlib.Path` instances.

    """
    return _glob(str(pathname), **kwargs)


def make_pdf(
    tmpdir=TEMP / "tex",
    publish=True,
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

    # Set up a temp directory and generate the stylesheet
    shutil.rmtree(tmpdir, ignore_errors=True)
    tmpdir.mkdir(parents=True)
    with open(tmpdir / f"{style}", "w") as f:
        sty = env.get_template(f"{style}.jinja").render(**jinja_kwargs)
        print(sty, file=f)

    # Copy LaTeX stylesheets & helper files to tmpdir
    for file in glob(ROOT / "tex" / "*"):
        shutil.copy(file, tmpdir)

    # Copy user files to tmpdir
    for file in glob(USER / "tex" / "*"):
        shutil.copy(file, tmpdir)

    # Debug mode?
    if settings.verbose:
        tectonic_args += ["--print"]

    # Build
    subprocess.check_output(
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
