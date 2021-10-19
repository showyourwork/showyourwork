rule pdf:
    """
    Build the article PDF.
    
    """
    message:
        "Building pdf..."
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
        report("ms.pdf", category="Article")
    params:
        verbose=config["verbose"],
        TEMP=relpaths.temp,
        TEX=relpaths.tex,
        SYWTEXFILE=files.tmp_syw,
        TECTONIC=config["tectonic_cmd"]
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/pdf.py"
