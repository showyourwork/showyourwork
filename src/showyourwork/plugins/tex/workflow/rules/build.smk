build_dir = SYW__WORK_PATHS / "build"

def _build_dependendencies_for(doc):
    deps_func = get_document_dependencies(doc)
    def impl(*_):
        deps = deps_func()
        return [build_dir / dep for dep in deps]
    return impl

for doc in SYW__DOCUMENTS:
    doc_dir = Path(doc).parent
    name = paths.path_to_rule_name(doc)
    pdf = SYW__WORK_PATHS.output / Path(doc).with_suffix(".pdf")

    rule:
        """
        Compile the document using ``tectonic``.
        """
        name:
            f"sywplug__tex_build_{name}"
        input:
            dependencies=_build_dependendencies_for(doc),
            document=build_dir / doc,
            style=build_dir / doc_dir / "showyourwork.tex",
            classfile=build_dir / doc_dir / "showyourwork.sty",
        output:
            pdf,
            output_directory=directory(SYW__WORK_PATHS.output),
        conda:
            SYWPLUG__TEX_RESOURCE("envs", "tectonic.yml")
        shell:
            "tectonic "
            "--chatter minimal "
            "--synctex "
            "--keep-logs --keep-intermediates "
            "--outdir {output.output_directory:q} "
            "{input.document:q} "

    rule:
        """
        Copy the output PDF from the output directory to the same directory as
        the source file.
        """
        name:
            f"sywplug__tex_output_build_{name}"
        input:
            pdf
        output:
            Path(doc).with_suffix(".pdf")
        run:
            utils.copy_file_or_directory(input[0], output[0])
