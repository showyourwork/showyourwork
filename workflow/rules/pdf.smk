rule pdf:
    input:
        preprocess((paths.temp / "config.json").relative_to(paths.user).as_posix())
    output:
        touch("ms.pdf")