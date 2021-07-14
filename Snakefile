configfile: "config.yml"


include: ".cortex/functions.py"


rule pdf:
    input: 
        "tex/ms.tex", 
        "tex/bib.bib",
        ".cortex/data/meta.json", 
        figures,
        test_results,
        "tests/tests.status"
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


rule foo:
    input:
        test_results
    output:
        "tests/tests.status"
    run:
        failures = 0
        for file in input:
            with open(file, "r") as f:
                if "ctxTestFailed" in f.read():
                    failures += 1
        if failures == 0:
            badge = r"\\\def\\\ctxTestsBadge{\\\color{ctxTestPassed}\\\faCheck}"
        else:
            badge = r"\\\def\\\ctxTestsBadge{\\\color{ctxTestFailed}\\\faTimes}"
        shell("echo {badge} > {output[0]}")


checkpoint script_meta:
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