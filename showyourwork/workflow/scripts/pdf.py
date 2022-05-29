"""
Compiles the article manuscript into a PDF.

"""
from showyourwork import paths
from showyourwork.tex import compile_tex
import sys
import shutil
from pathlib import Path
from jinja2 import Environment, BaseLoader


if __name__ == "__main__":

    # Snakemake config (available automagically)
    config = snakemake.config  # type:ignore

    # Metadata file jinja template
    TEMPLATE = r"""
    ((* if github_actions *))
    \OnGithubActionstrue
    ((* else *))
    \OnGithubActionsfalse
    ((* endif *))
    \def\syw@url{((- git_url -))}
    \def\syw@sha{((- git_sha -))}
    \def\syw@runid{((- github_runid -))}

    ((* for key, value in labels.items() *))
    \addvalue{((- key -))}{((- value -))}
    ((* endfor *))

    % Check if the Git tag is set (i.e. is it an empty string)
    ((* if sha_tag_header != "" *))
    \newcommand{\gitHeader}{((- sha_tag_header -))}
    ((* endif *))
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
    with open(str(Path(snakemake.config["stylesheet_meta_file"])), "w") as f:
        meta = ENV.from_string(TEMPLATE).render(**snakemake.config)
        print(meta, file=f)

    # Build the paper
    compile_tex(
        snakemake.config,
        output_dir=paths.user().compile,
        stylesheet=paths.showyourwork().resources / "styles" / "build.tex",
    )

    # Copy the PDF to the user dir
    shutil.copy(
        str(paths.user().compile / (snakemake.config["ms_name"] + ".pdf")),
        str(Path(snakemake.config["ms_pdf"])),
    )
