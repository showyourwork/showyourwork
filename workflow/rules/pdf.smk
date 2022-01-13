rule pdf:
    """
    Compile the manuscript into the article PDF.

    """
    input:
        config["ms_tex"],
        config["pdf_dependencies"],
        "showyourwork.yml"
    output:
        config["ms_pdf"],
        temp(config["tex_files_out"])
    conda:
        "../envs/main.yml"
    script:
        "../scripts/pdf.py"