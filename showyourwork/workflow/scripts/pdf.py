"""
Compiles the article manuscript into a PDF.

"""
import shutil
from pathlib import Path

from jinja2 import BaseLoader, Environment

from showyourwork import paths
from showyourwork.tex import compile_tex

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
    with open(str(Path(snakemake.config["stylesheet_meta_file"])), "w") as f:
        meta = ENV.from_string(TEMPLATE).render(**snakemake.config)
        print(meta, file=f)

    # Custom tectonic args for this step
    args = []
    if snakemake.config["synctex"]:
        args = ["--synctex"]

    # Build the paper
    compile_tex(
        snakemake.config,
        output_dir=paths.user().compile,
        stylesheet=paths.showyourwork().resources / "styles" / "build.tex",
        args=args,
    )

    # Copy the PDF to the user dir
    shutil.copy(
        str(paths.user().compile / (snakemake.config["ms_name"] + ".pdf")),
        str(Path(snakemake.config["ms_pdf"])),
    )

    # Copy the synctex file to the user dir
    if snakemake.config["synctex"]:
        path = paths.user().compile / (
            snakemake.config["ms_name"] + ".synctex.gz"
        )
        shutil.copy(str(path), snakemake.config["ms_name"] + ".synctex.gz")
