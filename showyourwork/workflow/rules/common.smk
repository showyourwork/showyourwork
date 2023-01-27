from showyourwork import paths

def temporary_tex_files():
    root = paths.user().compile
    files = [
        config["ms_tex"], config["stylesheet"], config["stylesheet_meta_file"]
    ] + list(config["tex_files_in"])
    return [
        (root / Path(f).name).as_posix()
        for f in files
    ]
