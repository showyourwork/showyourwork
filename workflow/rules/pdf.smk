rule pdf:
    input:
        f"src/tex/{config['ms_name']}.tex",
        "showyourwork.yml"
    output:
        f"{config['ms_name']}.pdf",
        temp(config["tex_files_out"])
    conda:
        "../envs/main.yml"
    script:
        "../scripts/pdf.py"