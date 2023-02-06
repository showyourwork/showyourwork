from showyourwork import paths

rule _syw_default_copy_manuscript:
    input:
        paths.repo(config).manuscript
    output:
        paths.work(config).manuscript
    shell:
        """
        cp "{input}" "{output}"
        """

repo_path = paths.repo(config).root
build_path = paths.work(config).build

for static_file in config.get("files", {}).get("static", []):
    rule:
        name:
            f"_syw_default_copy_static_{static_file}"
        input:
            repo_path / static_file
        output:
            directory(build_path / static_file) if (repo_path / static_file).is_dir() else build_path / static_file
        run:
            import shutil
            output[0].parent.mkdir(parents=True, exist_ok=True)            
            if input[0].is_dir():
                shutil.copytree(input[0], output[0])
            else:
                shutil.copyfile(input[0], output[0])


