from showyourwork import paths, utils

repo_path = paths.repo(config).root
work_path = paths.work(config).root

for doc in config.get("documents", ["ms.tex"]):
    doc_dir = Path(doc).parent
    name = Path(doc).name

    rule:
        name:
            f"syw__copy_doc_{doc}"
        input:
            repo_path / doc
        output:
            work_path / doc
        run:
            utils.copy_file_or_directory(input[0], output[0])
