"""
Defines the rule ``syw__ar5ivist``.

Runs ar5ivist to generate an HTML version of the manuscript.

"""
rule:
    """
    Generate a tarball for arXiv submission.

    """
    name:
        "syw__ar5ivist"
    message:
        "Generating the arXiv tarball..."
    input:
        temporary_tex_files(),
        compile_dir=paths.user().compile.as_posix()
    output:
        (paths.user().compile_html / "index.html").as_posix(),
    params:
        html_dir=paths.user().compile_html.as_posix(),
    shell:
        """
        set -e
        cd {input.compile_dir}
        mkdir -p syw_html
        docker run -v "$(pwd)":/docdir -w /docdir \
                    --user "$(id -u):$(id -g)" \
                    latexml/ar5ivist:2301.01 --source=ms.tex --destination=syw_html/index.html
        mv syw_html {params.html_dir}
        """
