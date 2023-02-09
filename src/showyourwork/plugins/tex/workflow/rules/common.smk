from functools import partial
from pathlib import Path

from showyourwork import paths

plugin_id = "showyourwork.plugins.tex"
resource = partial(paths.package_data, plugin_id, "workflow")
work_directory = paths.work(config).root
manuscript_name = config.get("manuscript", "ms.tex")
manuscript_directory = Path(manuscript_name).parent

def local_or_provided_style(*_):
    """Get the path to the showyourwork.sty file. We prefer to use the one
    provided by the user if it exists, but will provide our own if not.
    """
    path = paths.repo(config).manuscript.parent / "showyourwork.sty"
    if path.is_file():
        return path
    else:
        return resource("resources", "showyourwork.sty")

rule sywplug__tex_copy_style:
    input:
        resource("resources", "{runtype}.tex", check=False)
    output:
        work_directory / "{runtype}" / manuscript_directory / "showyourwork.tex"
    run:
        import shutil
        dst = Path(output[0])
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(input[0], dst)

rule sywplug__tex_copy_class:
    input:
        local_or_provided_style
    output:
        work_directory / "{runtype}" / manuscript_directory / "showyourwork.sty"
    run:
        import shutil
        dst = Path(output[0])
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(input[0], dst)
