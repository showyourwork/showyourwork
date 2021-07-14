configfile: "config.yml"


include: ".cortex/functions.py"


rule pdf:
    input: 
        "tex/ms.tex", 
        "tex/bib.bib",
        ".cortex/data/meta.json", 
        "tests/tests.status",
        figures,
        test_results
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
        figure=figure_wildcards
    params:
        script=figure_script_base_name,
        cache=figure_cache,
        other=figure_other
    run:
        run_figure(params)


rule test:
    input: 
        "{test}"
    output: 
        "{test}.status"
    wildcard_constraints:
        test=test_wildcards
    run: 
        run_test(input, output)


rule test_status:
    input:
        test_results
    output:
        "tests/tests.status"
    run:
        run_test_status(input, output)


checkpoint script_info:
    input: 
        "tex/ms.tex"
    output: 
        ".cortex/data/scripts.json"
    run: 
        get_script_metadata()
        get_metadata()


rule user_info:
    output: 
        ".cortex/data/user.json"
    run: 
        get_user_metadata()


rule metadata:
    input: 
        ".cortex/data/user.json", scripts
    output: 
        ".cortex/data/meta.json"
    run: 
        get_metadata()