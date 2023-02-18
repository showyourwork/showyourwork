import shutil
from functools import partial
from pathlib import Path

from showyourwork import paths, utils

plugin_id = "showyourwork.plugins.tex"
resource = partial(paths.package_data, plugin_id, "workflow")
work_directory = paths.work(config).root
manuscript_name = config.get("manuscript", "ms.tex")
manuscript_directory = Path(manuscript_name).parent

# def split_run_path(wildcards):
#     path = Path(wildcards.run_path)
#     return path.parts[0], Path(path.parts[1:])

def local_or_provided_style(document):
    """Get the path to the showyourwork.sty file. We prefer to use the one
    provided by the user if it exists, but will provide our own if not.
    """
    path = paths.repo(config).root / document
    path = path.parent / "showyourwork.sty"
    if path.is_file():
        return path
    else:
        return resource("resources", "showyourwork.sty")

# rule sywplug__tex_copy_style:
#     input:
#         lambda wildcards: resource(
#             "resources", f"{split_run_path(wildcards)[0]}.tex", check=False
#         )
#     output:
#         work_directory / "{run_path}" / "showyourwork.tex"
#     run:
#         dst = Path(output[0])
#         dst.parent.mkdir(parents=True, exist_ok=True)
#         shutil.copyfile(input[0], dst)

# rule sywplug__tex_copy_class:
#     input:
#         local_or_provided_style
#     output:
#         work_directory / "{run_path}" / "showyourwork.sty"
#     run:
#         import shutil
#         dst = Path(output[0])
#         dst.parent.mkdir(parents=True, exist_ok=True)
#         shutil.copyfile(input[0], dst)
