from showyourwork import paths

def temporary_tex_files(root=paths.user().compile):
    files = [
        config["ms_tex"], config["stylesheet"]
    ] + list(config["tex_files_in"])
    if root == paths.user().compile:
        files.append(config["stylesheet_meta_file"])
    return [
        (root / Path(f).name).as_posix()
        for f in files
    ]
