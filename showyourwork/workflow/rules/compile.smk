"""
Defines the rule ``syw__compile`` to build the article PDF.

"""
from showyourwork import paths


rule:
    """
    Compile the manuscript into the article PDF.

    """
    name:
        "syw__compile"
    message:
        "Generating the article PDF..."
    input:
        config["ms_tex"],
        config["dependencies"][config["ms_tex"]],
        "showyourwork.yml"
    output:
        config["ms_pdf"],
        temp(config["tex_files_out"]),
        temp(config["stylesheet"]),
        temp(config["stylesheet_meta_file"]),
        directory(paths.user().compile.as_posix())
    script:
        "../scripts/pdf.py"