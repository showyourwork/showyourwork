"""
Rules in this Snakefile are responsible for parsing the dependencies of any TeX
documents in the project and writing them in a standard JSON format. This works
by redefining TeX macros like ``\includegraphics`` and ``\input`` to write to a
log file instead of including the file. After running this hacked TeX file
through ``tectonic``, the log file is parsed to extract the dependency
structure. Of note, this procedure depends on the TeX documents including the
``showyourwork`` package, and the workflow will fail if it is not included.
"""

from showyourwork.plugins.tex.dependencies import parse_dependencies

deps_dir = SYW__WORK_PATHS / "dependencies"
style_paths = set()

# We only define these rules for documents explicitly listed in the config file
# because we otherwise end up with ambigious rules for other TeX files. So here
# we're looping over all the document paths.
for doc in SYW__DOCUMENTS:
    doc_dir = Path(doc).parent
    name = paths.path_to_rule_name(doc)
    xml =  deps_dir / doc_dir / f"{Path(doc).with_suffix('').name}.dependencies.xml"

    # If multiple documents live within the same directory, we only want to copy
    # the style files once.
    if str(doc_dir) not in style_paths:
        style_paths.add(str(doc_dir))

        rule:
            """
            Copy the appropriate ``showyourwork`` style file to the dependencies
            work directory. In this case, we're using the ``dependencies.tex``
            style, which redefines all the relevant macros to log their inputs.
            """
            name:
                f"sywplug__tex_deps_style_{name}"
            input:
                SYWPLUG__TEX_RESOURCE(
                    "resources", "dependencies.tex", check=False
                )
            output:
                deps_dir / doc_dir / "showyourwork.tex"
            run:
                utils.copy_file_or_directory(input[0], output[0])

        rule:
            """
            Copy the ``showyourwork`` class file to the dependencies work
            directory. If the project contains a ``showyourwork.sty`` file in
            the same directory as the document, we use that instead of the
            standard one provided by showyourwork, allowing users to customize
            behavior.
            """
            name:
                f"sywplug__tex_deps_class_{name}"
            input:
                sywplug__tex_local_or_provided_style(doc)
            output:
                deps_dir / doc_dir / "showyourwork.sty"
            run:
                utils.copy_file_or_directory(input[0], output[0])

    rule:
        """
        Copy the document from the parent work directory to the dependencies
        work directory.
        """
        name:
            f"sywplug__tex_deps_doc_{name}"
        input:
            SYW__WORK_PATHS.root / doc
        output:
            deps_dir / doc
        run:
            utils.copy_file_or_directory(input[0], output[0])

    rule:
        """
        Run tectonic to produce the XML log file tracking all data dependencies.
        """
        name:
            f"sywplug__tex_deps_xml_{name}"
        input:
            manuscript=deps_dir / doc,
            style=deps_dir / doc_dir / "showyourwork.tex",
            classfile=deps_dir / doc_dir / "showyourwork.sty",
        output:
            xml
        conda:
            SYWPLUG__TEX_RESOURCE("envs", "tectonic.yml")
        shell:
            """
            tectonic                 \\
                --chatter minimal    \\
                --keep-logs          \\
                --keep-intermediates \\
                {input.manuscript:q}
            """

    rule:
        """
        Parse the XML log file to produce a JSON file with the correctly
        formatted dependency structure.
        """
        name:
            f"sywplug__tex_deps_{name}"
        input:
            xml
        output:
            SYW__WORK_PATHS.dependencies_for(doc)
        run:
            base_path = (SYW__REPO_PATHS.root / doc).parent
            parse_dependencies(input[0], output[0], base_path)

rule sywplug__tex_dependencies:
    """
    A dummy rule that can be used to produce the dependencies for all TeX
    documents without knowing their specific names.
    """
    input:
        [SYW__WORK_PATHS.dependencies_for(doc) for doc in SYW__DOCUMENTS]
    output:
        touch(SYW__WORK_PATHS.flag("tex_dependencies"))
