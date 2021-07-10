import jinja2
import os
import json
from pathlib import Path
import subprocess
import glob


# Paths
CORTEX = Path(__file__).parents[0].absolute()


def generate_sty():
    """
    Generate a new `styles/cortex.sty` file from a template.

    """
    # Special delimiters for LaTeX
    env = jinja2.Environment(
        block_start_string="((*",
        block_end_string="*))",
        variable_start_string="((-",
        variable_end_string="-))",
        comment_start_string="((=",
        comment_end_string="=))",
        trim_blocks=True,
        autoescape=False,
        loader=jinja2.FileSystemLoader(CORTEX),
    )

    # Load the metadata
    with open(CORTEX / "data" / "meta.json", "r") as f:
        meta = json.load(f)

    # Get the current git hash
    try:
        meta["repo"]["sha"] = (
            subprocess.check_output(["git", "rev-parse", "HEAD"])
            .decode()
            .replace("\n", "")
        )
    except:
        meta["repo"]["sha"] = meta["repo"]["branch"]

    # Generate the style file
    with open(CORTEX / "styles" / "cortex.sty", "w") as f:
        print(
            env.get_template(
                str(Path("templates") / "cortex.sty.jinja")
            ).render(**meta),
            file=f,
        )


if __name__ == "__main__":
    generate_sty()
