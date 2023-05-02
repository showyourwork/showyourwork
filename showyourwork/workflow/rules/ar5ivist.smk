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
        paths.user().html / "index.html",
    shell:
        """
        cd {input.compile_dir}
        docker run -v "$(pwd)":/docdir -w /docdir \
                    --user "$(id -u):$(id -g)" \
                    latexml/ar5ivist:2301.01 --source=ms.tex --destination=tmphtml/index.html
        mv tmphtml/* {output}/
        """
