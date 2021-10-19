localrules: class_file


rule class_file:
    """
    Copy the LaTeX class file to the ``src`` directory.

    """
    message:
        "Generating auxiliary class file `{output}`..."
    input:
        input_class_file,
    output:
        temp(posix(relpaths.tex / "{file}")),
    wildcard_constraints:
        file="|".join(
            [file for files in files.cls.values() for file in files]
        ),
    shell:
        "cp {input[0]} {output}"