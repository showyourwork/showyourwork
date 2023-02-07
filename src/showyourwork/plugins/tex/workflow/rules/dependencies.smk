from functools import partial

from showyourwork import paths
from showyourwork.plugins.tex.dependencies import parse_dependencies

plugin_id = "showyourwork.plugins.tex"
xml_directory = paths.work(config).plugin(plugin_id, "xml")
resource = partial(paths.package_data, plugin_id, "workflow")

#
# Rules for generating the XML dependency tree from a TeX manuscript
#
rule syw__plug_tex_copy_ms2xml:
    input:
        paths.work(config).manuscript
    output:
        xml_directory / "manuscript.tex"
    shell:
        """
        cp "{input}" "{output}"
        """

rule syw__plug_tex_copy_style2xml:
    input:
        resource("resources", "xmlstyle.tex")
    output:
        xml_directory / "showyourwork.tex"
    shell:
        """
        cp "{input}" "{output}"
        """

def local_or_provided_style(*args):
    """Get the path to the showyourwork.sty file. We prefer to use the one
    provided by the user if it exists, but will provide our own if not.
    """
    path = paths.repo(config).manuscript.parent / "showyourwork.sty"
    if path.is_file():
        return path
    else:
        return resource("resources", "showyourwork.sty")

rule syw__plug_tex_copy_class2xml:
    input:
        local_or_provided_style
    output:
        xml_directory / "showyourwork.sty"
    shell:
        """
        cp "{input}" "{output}"
        """

rule syw__plug_tex_xml_tree:
    input:
        manuscript=xml_directory / "manuscript.tex",
        style=xml_directory / "showyourwork.tex",
        classfile=xml_directory / "showyourwork.sty",
    output:
        xml_directory / "showyourwork.xml"
    conda: 
        resource("envs", "tectonic.yml")
    shell:
        """
        tectonic                  \\
            --chatter minimal     \\
            --keep-logs           \\
            --keep-intermediates  \\
            "{input.manuscript}"
        """

rule syw__plug_tex_dependencies:
    input:
        xml_directory / "showyourwork.xml"
    output:
        paths.work(config).dependencies
    run:
        parse_dependencies(input[0], output[0], paths.repo(config).manuscript.parent) 
