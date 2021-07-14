configfile: "config.yml"
debug = config["debug"]
quiet = config["quiet"]
ext = "|".join(config["figure_extensions"])


include: ".cortex/functions.py"


rule pdf:
    input: 
        "tex/ms.tex", 
        "tex/bib.bib",
        ".cortex/data/meta.json", 
        figures
    output: 
        "ms.pdf"
    run: 
        gen_pdf()


rule figure:
    input: 
        figure_script
    output: 
        "{figure}"
    wildcard_constraints:
        figure="figures/(.*?)\.{}".format(ext)
    params:
        script=figure_script_base_name,
        cache=figure_cache,
        other=figure_other
    run:
        run_figure(params)


checkpoint check_scripts:
    input: 
        "tex/ms.tex"
    output: 
        ".cortex/data/scripts.json"
    run: 
        get_script_metadata()


rule user_meta:
    output: 
        ".cortex/data/user.json"
    run: 
        get_user_metadata()


rule meta:
    input: 
        ".cortex/data/user.json", scripts
    output: 
        ".cortex/data/meta.json"
    run: 
        get_metadata()