sywplug_tex__rule_name = partial(
    utils.rule_name, plugin="showyourwork.plugins.tex"
)
sywplug_tex__resource = partial(
    package_data, "showyourwork.plugins.tex", "workflow"
)

def sywplug_tex__local_or_provided_style(document):
    """Get the path to the showyourwork.sty file. We prefer to use the one
    provided by the user if it exists, but will provide our own if not.
    """
    path = (SYW__REPO_PATHS.root / document).parent / "showyourwork.sty"
    if path.is_file():
        return path
    else:
        return sywplug_tex__resource("resources", "showyourwork.sty")


for base_path in [SYW__WORK_PATHS / "dependencies", SYW__WORK_PATHS / "build"]:
    slug = base_path.name

    rule:
        """
        Copy explicit dependencies to the working directory.
        """
        name:
            sywplug_tex__rule_name("copy", slug)
        input:
            "{file}"
        output:
            base_path / "{file}"
        run:
            utils.copy_file_or_directory(input[0], output[0])

    style_paths = set()
    for doc in SYW__DOCUMENTS:
        doc_dir = Path(doc).parent

        # If multiple documents live within the same directory, we only want to copy
        # the style files once.
        if str(doc_dir) not in style_paths:
            style_paths.add(str(doc_dir))

            rule:
                """
                Copy the appropriate ``showyourwork`` style file to the work
                directory.
                """
                name:
                    sywplug_tex__rule_name("style", slug, document=doc)
                input:
                    sywplug_tex__resource("resources", f"{slug}.tex")
                output:
                    base_path / doc_dir / "showyourwork.tex"
                run:
                    utils.copy_file_or_directory(input[0], output[0])

            rule:
                """
                Copy the ``showyourwork`` class file to the work directory. If
                the project contains a ``showyourwork.sty`` file in the same
                directory as the document, we use that instead of the standard
                one provided by showyourwork, allowing users to customize
                behavior.
                """
                name:
                    sywplug_tex__rule_name("class", slug, document=doc)
                input:
                    sywplug_tex__local_or_provided_style(doc)
                output:
                    base_path / doc_dir / "showyourwork.sty"
                run:
                    utils.copy_file_or_directory(input[0], output[0])

        rule:
            """
            Copy the document from the parent work directory to the dependencies
            work directory.
            """
            name:
                sywplug_tex__rule_name("doc", slug, document=doc)
            input:
                SYW__WORK_PATHS.root / doc
            output:
                base_path / doc
            run:
                utils.copy_file_or_directory(input[0], output[0])
