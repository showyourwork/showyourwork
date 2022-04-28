"""
Defines the rule ``syw__preprocess`` to parse the config and build the
workflow graph.

"""
from showyourwork import paths


rule:
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
    name:
        "syw__preprocess"
    message:
        "Setting up the workflow..."
    input:
        config["ms_tex"],
        "showyourwork.yml"
    output:
        config["config_json"],
        temp(config["tex_files_out"]),
        temp(config["stylesheet"]),
        directory(paths.user().preprocess.as_posix())
    script:
        "../scripts/preprocess.py"