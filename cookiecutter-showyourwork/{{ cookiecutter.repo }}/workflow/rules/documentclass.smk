def input_class_file(wildcards):
    """

    """
    checkpoints.class_name.get(**wildcards)
    with open(TEMP / "class_name", "r") as f:
        folder = f.read().replace("\n", "")
    return posix(WORKFLOW / "resources" / "classes" / folder / wildcards.file)


def class_files(wildcards):
    """

    """
    checkpoints.class_name.get(**wildcards)
    with open(TEMP / "class_name", "r") as f:
        folder = f.read().replace("\n", "")
    return [posix(TEX / file.name) for file in (WORKFLOW / "resources" / "classes" / folder).glob("*.*")]


checkpoint class_name:
    message:
        "Inferring document class..."
    output:
        temp(posix(TEMP / "class_name"))
    priority:
        100
    run:
        with open(TEX / "ms.tex", "r") as f:
            lines = f.readlines()
            for line in lines:
                match = re.match("[ \t]*\\\documentclass\[?.*?\]?\{(.*?)\}", line)
                if hasattr(match, "groups"):
                    name = match.groups()[0]
                    break
            else:
                raise ValueError("Unable to determine document class in `tex/ms.tex`.")
            if not TEMP.exists():
                os.mkdir(TEMP)
            with open(TEMP / "class_name", "w") as f:
                print(name, file=f)


rule class_file:
    message:
        "Generating auxiliary class file `{output}`..."
    input:
        input_class_file
    output:
        temp(posix(TEX / "{file}"))
    wildcard_constraints:
        file="|".join([file.name for file in (WORKFLOW / "resources" / "classes").glob("*/*")])
    shell:
        "cp {input[0]} {output}"
