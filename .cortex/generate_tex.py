import jinja2
import os
import json
from pathlib import Path


# Paths
CORTEX = Path(__file__).parents[0].absolute()
TEX = Path(__file__).parents[1].absolute() / "tex"


def generate_tex():
    """
    Generate a new `/tex/ms.tex` file from a template.

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

    # Generate the tex file
    with open(TEX / "ms.tex", "w") as f:
        print(
            env.get_template(str(Path("templates") / "ms.tex.jinja")).render(
                **meta
            ),
            file=f,
        )


if __name__ == "__main__":
    generate_tex()
