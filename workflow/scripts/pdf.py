import sys
import shutil


# Import utils
sys.path.insert(1, snakemake.config["workflow_abspath"])
from utils import paths, compile_tex

compile_tex(
    args=[
        "--chatter",
        "minimal",
        "--keep-logs",
        "--keep-intermediates",
        "-o",
        str(paths.build),
    ],
    stylesheet=paths.resources / "styles" / "build.tex",
    config=snakemake.config,
)

shutil.copy(
    str(paths.build / (snakemake.config["ms_name"] + ".pdf")),
    str(snakemake.config["ms_pdf"]),
)
