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
        using_preprocessing_script = (
            config["preprocess_arxiv_script"] is not None
        )
        if using_preprocessing_script:
            # First, copy to a temporary directory, excluding files
            shutil.copytree(paths.user().compile, tmpdir, dirs_exist_ok=True)
            # Run the preprocessing script, if provided:
            script = config["preprocess_arxiv_script"]
            subprocess.run([script, tmpdir], check=True)

        src_dir = (
            Path(tmpdir)
            if using_preprocessing_script
            else paths.user().compile
        )
        # Tar up everything in the src/tex directory
        with tarfile.open("arxiv.tar.gz", "w:gz") as tarball:
            for file in src_dir.rglob("*"):
                if file.name not in exclude and file.suffix != ".bib":
                    tarball.add(file, arcname=file.name)
