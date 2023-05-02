"""
Wraps the :doc:`pdf` script to build the article PDF, then
tars everything in the ``src/tex`` directory into the tarball
``arxiv.tar.gz`` for easy arXiv submission.

"""
import shutil
import subprocess
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

from showyourwork import paths

if __name__ == "__main__":
    # Snakemake config (available automagically)
    config = snakemake.config  # type:ignore

    # File names to exclude
    ms_name = config["ms_name"]
    exclude = [
        ".gitignore",
        f"{ms_name}.pdf",
        f"{ms_name}.aux",
        f"{ms_name}.blg",
        f"{ms_name}.log",
        f"{ms_name}.out",
        f"{ms_name}.fls",
        f"{ms_name}.synctex.gz",
        f"{ms_name}.fdb_latexmk",
    ]

    with TemporaryDirectory() as tmpdir:
        # First, copy to a temporary directory, excluding files
        shutil.copytree(paths.user().compile, tmpdir, dirs_exist_ok=True)

        if config["preprocess_arxiv"]["enabled"]:
            # Run the preprocessing script, if provided:
            script = config["preprocess_arxiv"]["script"]
            if not script:
                raise ValueError(
                    "Preprocessing script not provided in config "
                    "(`preprocess_arxiv: script: '...'`)."
                )
            subprocess.run([script, tmpdir], check=True)

        # Tar up everything in the src/tex directory
        with tarfile.open("arxiv.tar.gz", "w:gz") as tarball:
            for file in Path(tmpdir).rglob("*"):
                if file.name not in exclude:
                    tarball.add(file, arcname=file.name)
