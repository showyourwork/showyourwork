rule arxiv:
    message:
        "Building arxiv tarball..."
    input:
        POSIX(TEX / "{}.tex".format(SYWTEXFILE)),
        [POSIX(file) for file in TEX.glob("*.bib")],
        POSIX(TEMP / "meta.json"),
        POSIX(TEX / "showyourwork.sty"),
        AUXFILES,
        class_files,
        figures,
    output:
        "arxiv.tar.gz"
    params:
        verbose=verbose,
        figexts=figexts,
        arxiv_tarball_exclude=arxiv_tarball_exclude,
        TEMP=TEMP,
        FIGURES=FIGURES,
        TEX=TEX,
        SYWTEXFILE=SYWTEXFILE,
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/arxiv.py"
