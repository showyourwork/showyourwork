"""
Compiles the article manuscript into a PDF.

"""

import shutil
from pathlib import Path

from jinja2 import BaseLoader, Environment

from showyourwork import paths


if __name__ == "__main__":
    # Snakemake config (available automagically)
    config = snakemake.config  # type:ignore

    # Create the compile directory
    compile_dir = Path(snakemake.output.compile_dir)

    # Copy over the source files
    shutil.copytree(paths.user().tex, compile_dir, dirs_exist_ok=True)

    if snakemake.params.metadata:
        # Metadata file jinja template
        TEMPLATE = r"""
        ((* if github_actions *))
        \OnGithubActionstrue
        ((* else *))
        \OnGithubActionsfalse
        ((* endif *))
        \def\syw@url{((- git_url -))}
        \def\syw@slug{((- git_slug -))}
        \def\syw@sha{((- git_sha -))}
        \def\syw@runid{((- github_runid -))}

        ((* if stamp.enabled *))
        \AddStamptrue
        ((* else *))
        \AddStampfalse
        ((* endif *))
        \newcommand{\syw@stampX}{((- stamp.xpos -))}
        \newcommand{\syw@stampY}{((- stamp.ypos -))}
        \newcommand{\syw@stampSize}{((- stamp.size -))}
        \newcommand{\syw@stampAngle}{((- stamp.angle -))}
        \newcommand{\syw@stampText}{((- stamp.text -))}
        \newcommand{\syw@stampVersion}{((- stamp.version -))}

        ((* for key, value in labels.items() *))
        \addvalue{((- key -))}{((- value -))}
        ((* endfor *))

        ((* for key, value in variables.items() *))
        \addvalue{((- key -))}{((- value -))}
        ((* endfor *))
        """

        # Custom jinja environment for LaTeX
        ENV = Environment(
            block_start_string="((*",
            block_end_string="*))",
            variable_start_string="((-",
            variable_end_string="-))",
            comment_start_string="((=",
            comment_end_string="=))",
            trim_blocks=True,
            autoescape=False,
            loader=BaseLoader(),
        )

        # Generate the stylesheet metadata file
        p = compile_dir / Path(config["stylesheet_meta_file"]).name
        with open(p, "w") as f:
            meta = ENV.from_string(TEMPLATE).render(**config)
            print(meta, file=f)

    # Copy over stylesheet
    p = compile_dir / Path(config["stylesheet"]).name
    shutil.copy(snakemake.input.stylesheet, p)

    # Copy over TeX auxiliaries
    for file in config["tex_files_in"]:
        dst = compile_dir / Path(file).name
        if not dst.exists():
            shutil.copy(file, dst)
