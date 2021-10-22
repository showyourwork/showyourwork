rule arxiv:
    """
    Builds a tarball of the article PDF and all output for posting to the arXiv.

    """
    message:
        "Building arxiv tarball..."
    input:
        posix(relpaths.tex / "{}.tex".format(files.tmp_syw)),
        [posix(file) for file in relpaths.tex.glob("*.bib")],
        posix(relpaths.temp / "meta.json"),
        posix(relpaths.tex / "showyourwork.sty"),
        files.aux,
        files.tectonic,
        class_files,
        figures,
    output:
        "arxiv.tar.gz"
    params:
        verbose=config["verbose"],
        figexts=config["figexts"],
        arxiv_tarball_exclude=config["arxiv_tarball_exclude"],
        ZENODO_FILES=files.zenodo_files_manual + files.zenodo_files_auto,
        TEMP=relpaths.temp,
        FIGURES=relpaths.figures,
        TEX=relpaths.tex,
        SYWTEXFILE=files.tmp_syw,
        TECTONIC=config["tectonic_cmd"]
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/arxiv.py"
