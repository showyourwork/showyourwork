rule sywtexfile:
    message:
        "Writing temporary tex file..."
    input:
        posix(TEX / "ms.tex")
    output:
        temp(posix(TEX / "{}.tex".format(SYWTEXFILE)))
    run:
        with open(TEX / "ms.tex", "r") as f:
            lines = f.readlines()
        for idx, line in enumerate(lines):
            if line.startswith(r"\documentclass"):
                lines = (
                    lines[: idx + 1]
                    + [r"\usepackage{showyourwork}" + "\n"]
                    + lines[idx + 1 :]
                )
                break
        else:
            raise ValueError(r"Missing `\documentclass` in file `tex/ms.tex`.")
        with open(TEX / "{}.tex".format(SYWTEXFILE), "w") as f:
            f.writelines(lines)


rule stylesheet:
    input:
        posix(TEMP / "meta.json")
    output:
        temp(posix(TEX / "showyourwork.sty"))
    run:
        env = jinja2.Environment(
            block_start_string="((*",
            block_end_string="*))",
            variable_start_string="((-",
            variable_end_string="-))",
            comment_start_string="((=",
            comment_end_string="=))",
            trim_blocks=True,
            autoescape=False,
            loader=jinja2.FileSystemLoader(WORKFLOW / "resources" / "templates"),
        )
        with open(TEMP / "meta.json", "r") as f:
            jinja_kwargs = json.load(f)
        with open(TEX / "showyourwork.sty", "w") as f:
            sty = env.get_template("showyourwork.sty").render(**jinja_kwargs)
            print(sty, file=f)


rule pdf:
    message:
        "Building pdf..."
    input:
        posix(TEX / "{}.tex".format(SYWTEXFILE)),
        [posix(file) for file in TEX.glob("*.bib")],
        TEMP / "meta.json",
        posix(TEX / "showyourwork.sty"),
        class_files,
        aux_files,
        figures
    output:
        temp(posix(TEMP / "{}.pdf".format(SYWTEXFILE))),
        "ms.pdf"
    run:
        tectonic_args = ["-o", TEMP]
        if verbose:
            tectonic_args += ["--print"]
        else:
            tectonic_args += ["--chatter", "minimal"]
        subprocess.check_call(
            ["tectonic"] + tectonic_args + [TEX / "{}.tex".format(SYWTEXFILE)]
        )
        shutil.copy(TEMP / "{}.pdf".format(SYWTEXFILE), "ms.pdf")
