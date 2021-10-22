localrules: xml


rule xml:
    """
    Generates the article XML tree. Specifically, builds the article
    using ``tectonic``, but re-defines ``figure``, ``caption``, and ``label``
    commands to print XML tags to a special log file. This way, we can
    use LaTeX to construct a full XML tree of the document for us, without
    any need for parsing the TeX file ourselves.
    This XML tree is then used to determine relationships between the figure
    scripts and the figure files.
    
    """
    message:
        "Generating XML tree..."
    input:
        class_files,
        files.aux,
        posix(relpaths.tex / "sywxml.sty"),
        posix(relpaths.tex / "{}.tex".format(files.tmp_xml)),
        files.tectonic
    output:
        posix(relpaths.temp / "showyourwork.xml"),
        temp(posix(relpaths.temp / "{}.pdf".format(files.tmp_xml))),
    params:
        verbose=config["verbose"],
        TEMP=relpaths.temp,
        TEX=relpaths.tex,
        TMPTEXFILE=files.tmp_xml,
        TECTONIC=config["tectonic_cmd"]
    conda:
        "../envs/environment.yml"
    script:
        "../scripts/xml.py"