from showyourwork import paths
from showyourwork.plugins.tex.dependencies import parse_dependencies

plugin_id = "showyourwork.plugins.tex"
repo_directory = paths.repo(config).root
work_directory = paths.work(config).root
build_directory = paths.work(config).build
output_directory = paths.work(config).output
resource = partial(paths.package_data, plugin_id, "workflow")

def _repo_to_build(path):
    return build_directory / Path(path).relative_to(repo_directory)

def _work_to_build(path):
    return build_directory / Path(path).relative_to(work_directory)

def _repo_paths(*_):
    return get_manuscript_dependencies()

def _work_paths(*_):
    return [manuscript_path]

def _build_dependendencies(*_):
    deps = get_manuscript_dependencies()
    return [_repo_to_build(dep) for dep in deps]

# rule sywplug__tex_copy_manuscript_dependencies:
#     input:
#         rules.syw__dag.output,
#         ensure_manuscript_dependencies,
#         paths=_repo_paths,
#         manuscript=paths.work(config).manuscript,
#     output:
#         directory(build_directory),
#         manuscript=build_directory / config.get("manuscript", "ms.tex"),
#     run:
#         import shutil

#         for src in input.paths:
#             dst = _repo_to_build(src)
#             dst.parent.mkdir(parents=True, exist_ok=True)
#             if Path(src).is_file():
#                 shutil.copy(src, dst)
#             else:
#                 shutil.copytree(src, dst)

#         src = input.manuscript
#         dst = Path(output.manuscript)
#         dst.parent.mkdir(parents=True, exist_ok=True)
#         shutil.copy(src, dst)

rule sywplug__tex_copy_files_to_build:
    input:
        rules.syw__dag.output,
        ensure_manuscript_dependencies,
        source=paths.repo(config).root / "{file}"
    output:
        build_directory / "{file}"
    run:
        import shutil
        dst = Path(output[0])
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(input.source, dst)

rule sywplug__tex_build:
    input:
        dependencies=_build_dependendencies,
        manuscript=build_directory / manuscript_name,
        style=build_directory / manuscript_directory / "showyourwork.tex",
        classfile=build_directory / manuscript_directory / "showyourwork.sty",
    output:
        output_directory / Path(config.get("manuscript", "ms.tex")).with_suffix(".pdf").name,
        output_directory=directory(output_directory),
    conda: 
        resource("envs", "tectonic.yml")
    shell:
        """
        tectonic                                 \\
            --chatter minimal                    \\
            --keep-logs                          \\
            --keep-intermediates                 \\
            --synctex                            \\
            --outdir {output.output_directory:q} \\
            {input.manuscript:q}
        """
