rule arxiv_tarball:
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
        tar=True,
        TEMP=TEMP,
        FIGURES=FIGURES,
        TEX=TEX,
        SYWTEXFILE=SYWTEXFILE,
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/arxiv.py"

rule arxiv_folder:
    message:
        "Building arxiv folder..."
    input:
        POSIX(TEX / "{}.tex".format(SYWTEXFILE)),
        [POSIX(file) for file in TEX.glob("*.bib")],
        POSIX(TEMP / "meta.json"),
        POSIX(TEX / "showyourwork.sty"),
        AUXFILES,
        class_files,
        figures,
    output:
        directory("arxiv")
    params:
        verbose=verbose,
        figexts=figexts,
        tar=False,
        TEMP=TEMP,
        FIGURES=FIGURES,
        TEX=TEX,
        SYWTEXFILE=SYWTEXFILE,
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/arxiv.py"
