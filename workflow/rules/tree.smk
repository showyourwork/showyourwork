rule tree:
    """
    Generate the article dependency tree. 
    
    This rule builds the article using ``tectonic``, but re-defines ``figure``, 
    ``caption``, and ``label`` commands to print XML tags to a special log file. 
    This way, we can use TeX to construct a full XML tree of the document for us, 
    without any need for parsing the TeX file ourselves. This XML tree is then 
    used to determine relationships between the figure scripts and the figure 
    files.
    
    """
    input:
        f"src/tex/{config['ms_name']}.tex"
    output:
        (paths.temp / "tree.json").as_posix(),
        temp((paths.temp / "showyourwork.xml").as_posix()),
        temp((paths.temp / f"{config['ms_name']}.pdf").as_posix()),
        [temp(f"src/tex/{Path(file).name}") for file in config["tex_files"]]
    conda:
        "../envs/setup.yml"
    script:
        "../scripts/tree.py"