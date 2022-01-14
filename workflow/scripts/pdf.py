import sys
import shutil
from pathlib import Path
from jinja2 import Environment, BaseLoader


# Import utils
sys.path.insert(1, snakemake.config["workflow_abspath"])
from utils import paths, compile_tex


# Metadata file jinja template
TEMPLATE = r"""
((* if github_actions *))
\OnGithubActionstrue
((* else *))
\OnGithubActionsfalse
((* endif *))
\def\syw@url{((- git_url -))}
\def\syw@sha{((- git_sha-))}
\def\syw@runid{((- github_runid-))}

((* for key, value in labels.items() *))
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


# Build the paper
compile_tex(
    args=[
        "--chatter",
        "minimal",
        "--keep-logs",
        "--keep-intermediates",
        "-o",
        str(paths.compile),
    ],
    stylesheet=paths.resources / "styles" / "build.tex",
    config=snakemake.config,
)


# Copy the PDF to the user dir
shutil.copy(
    str(paths.compile / (snakemake.config["ms_name"] + ".pdf")),
    str(Path(snakemake.config["ms_pdf"])),
)
