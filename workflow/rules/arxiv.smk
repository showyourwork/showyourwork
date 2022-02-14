from utils import paths


rule:
    """
    Generate a tarball for arXiv submission.

    """
    name:
        "syw__arxiv"
    message:
        "Generating the arXiv tarball..."
    input:
        config["ms_tex"],
        config["dependencies"][config["ms_tex"]],
        "showyourwork.yml"
    output:
        "arxiv.tar.gz",
        config["ms_pdf"],
        temp(config["tex_files_out"]),
        temp(config["stylesheet"]),
        temp(config["stylesheet_meta_file"]),
        directory(paths.compile.as_posix())
    conda:
        "../envs/main.yml"
    script:
        "../scripts/arxiv.py"