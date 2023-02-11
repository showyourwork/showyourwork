from showyourwork import paths
from showyourwork.plugins.tex.dependencies import parse_dependencies

dependencies_directory = paths.work(config) / "dependencies"

rule sywplug__tex_copy_manuscript_for_dependencies:
    input:
        paths.work(config).manuscript
    output:
        dependencies_directory / manuscript_name
    run:
        import shutil
        dst = Path(output[0])
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(input[0], dst)

rule sywplug__tex_xml_tree:
    input:
        manuscript=dependencies_directory / config.get("manuscript", "ms.tex"),
        style=dependencies_directory / manuscript_directory / "showyourwork.tex",
        classfile=dependencies_directory / manuscript_directory / "showyourwork.sty",
    output:
        dependencies_directory / manuscript_directory / "showyourwork.xml"
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

rule sywplug__tex_dependencies:
    input:
        rules.sywplug__tex_xml_tree.output[0]
    output:
        paths.work(config).dependencies
    run:
        parse_dependencies(input[0], output[0], paths.repo(config).manuscript.parent)
