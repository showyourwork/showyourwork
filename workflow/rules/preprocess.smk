from utils import paths


rule preprocess:
    """
    Generate a `config.json` file for the main build.
    
    This rule builds the article using ``tectonic``, but re-defines ``figure``, 
    ``caption``, and ``label`` commands to print XML tags to a special log file. 
    This way, we can use TeX to construct a full XML tree of the document for us, 
    without any need for parsing the TeX file ourselves. This XML tree is then 
    used to determine relationships between the figure scripts and the figure 
    files.

    This rule also assembles information about the datasets and other script
    dependencies, as well as metadata about the git repo. It then packages
    all this up alongside the user's config settings into the file
    `config.json`, which is used as input to the main `showyourwork`
    workflow.
    
    """
    message:
        "Setting up the workflow..."
    input:
        config["ms_tex"],
        "showyourwork.yml"
    output:
        config["config_json"],
        temp(config["tex_files_out"]),
        temp(config["stylesheet"]),
        temp(directory(paths.preprocess.as_posix()))
    conda:
        "../envs/main.yml"
    script:
        "../scripts/preprocess.py"