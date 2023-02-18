from pathlib import Path
from showyourwork import paths, utils
from showyourwork.plugins.tex.dependencies import parse_dependencies

work_paths = paths.work(config)
deps_dir = work_paths / "dependencies"
style_paths = set()
documents = config.get("documents", ["ms.tex"])

for doc in documents:

    doc_dir = Path(doc).parent
    name = Path(doc).name
    xml =  deps_dir / doc_dir / f"{Path(doc).stem}.dependencies.xml"

    if str(doc_dir) not in style_paths:
        style_paths.add(str(doc_dir))

        rule:
            name:
                f"sywplug__tex_deps_style_{name}"
            input:
                resource("resources", "dependencies.tex", check=False)
            output:
                deps_dir / doc_dir / "showyourwork.tex"
            run:
                utils.copy_file_or_directory(input[0], output[0])

        rule:
            name:
                f"sywplug__tex_deps_class_{name}"
            input:
                local_or_provided_style(doc)
            output:
                deps_dir / doc_dir / "showyourwork.sty"
            run:
                utils.copy_file_or_directory(input[0], output[0])

    rule:
        name:
            f"sywplug__tex_deps_doc_{name}"
        input:
            work_paths.root / doc
        output:
            deps_dir / doc
        run:
            utils.copy_file_or_directory(input[0], output[0])

    rule:
        name:
            f"sywplug__tex_deps_xml_{name}"
        input:
            manuscript=deps_dir / doc,
            style=deps_dir / doc_dir / "showyourwork.tex",
            classfile=deps_dir / doc_dir / "showyourwork.sty",
        output:
            xml
        conda:
            resource("envs", "tectonic.yml")
        shell:
            """
            tectonic                 \\
                --chatter minimal    \\
                --keep-logs          \\
                --keep-intermediates \\
                {input.manuscript:q}
            """

    rule:
        name:
            f"sywplug__tex_deps_{name}"
        input:
            xml
        output:
            work_paths.dependencies_for(doc)
        run:
            base_path = (paths.repo(config).root / doc).parent
            parse_dependencies(input[0], output[0], base_path)

rule sywplug__tex_dependencies:
    input:
        [work_paths.dependencies_for(doc) for doc in documents]
    output:
        touch(work_paths.flag("tex_dependencies"))
