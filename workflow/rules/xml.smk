localrules: xml


rule xml:
    """
    Generate the article XML tree.
    
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