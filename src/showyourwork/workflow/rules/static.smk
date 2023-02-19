from showyourwork import paths, utils

repo_path = paths.repo(config).root
build_path = paths.work(config).build

for static_file in config.get("static", []):
    rule:
        name:
            f"syw__static_{paths.path_to_rule_name(static_file)}"
        input:
            repo_path / static_file
        output:
            directory(build_path / static_file) if (repo_path / static_file).is_dir() else build_path / static_file
        run:
            utils.copy_file_or_directory(input[0], output[0])
