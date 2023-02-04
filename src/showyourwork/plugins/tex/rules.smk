from showyourwork import paths
from showyourwork.plugins.tex import add_to_preamble

plugin_id = "showyourwork.plugins.tex"
xml_directory = paths.work(config).plugin(plugin_id, "xml")

#
# Rules for generating the XML dependency tree from a TeX manuscript
#
rule _syw_plug_tex_copy_ms2xml:
    input:
        paths.repo(config).manuscript
    output:
        xml_directory / "manuscript.tex"
    shell:
        """
        cp "{input}" "{output}"
        """

rule _syw_plug_tex_copy_style2xml:
    input:
        paths.package_data("resources", "xmlstyle.tex")
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
        return paths.package_data("resources", "showyourwork.sty")

rule _syw_plug_tex_copy_class2xml:
    input:
        local_or_provided_style
    output:
        xml_directory / "showyourwork.sty"
    shell:
        """
        cp "{input}" "{output}"
        """

rule _syw_plug_tex_xml_tree:
    input:
        manuscript=xml_directory / "manuscript.tex",
        style=xml_directory / "showyourwork.tex",
        classfile=xml_directory / "showyourwork.sty",
    output:
        xml_directory / "showyourwork.xml"
    conda: 
        paths.package_data(plugin_id, "resources", "tectonic.yml")
    shell:
        """
        tectonic                  \\
            --chatter minimal     \\
            --keep-logs           \\
            --keep-intermediates  \\
            "{input.manuscript}"
        """
