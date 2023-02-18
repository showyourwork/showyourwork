"""
Wraps the :doc:`pdf` script to build the article PDF, then
tars everything in the ``src/tex`` directory into the tarball
``arxiv.tar.gz`` for easy arXiv submission.

"""
import tarfile

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

    # Tar up everything in the src/tex directory
    with tarfile.open("arxiv.tar.gz", "w:gz") as tarball:
        for file in paths.user().compile.rglob("*"):
            if file.name not in exclude:
                tarball.add(file, arcname=file.name)
