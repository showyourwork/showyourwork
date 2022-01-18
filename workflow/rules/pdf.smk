rule pdf:
    """
    Compile the manuscript into the article PDF.

    """
    message:
        "Generating the article PDF..."
    input:
        config["ms_tex"],
        config["pdf_dependencies"],
        "showyourwork.yml"
    output:
        config["ms_pdf"],
        temp(config["tex_files_out"]),
        temp(config["stylesheet"]),
        temp(config["stylesheet_meta_file"]),
        temp(directory(paths.compile.as_posix()))
    conda:
        "../envs/main.yml"
    script:
        "../scripts/pdf.py"