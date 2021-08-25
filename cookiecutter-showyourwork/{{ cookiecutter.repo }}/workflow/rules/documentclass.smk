import os
import re


def input_class_file(wildcards):
    """

    """
    checkpoints.class_name.get(**wildcards)
    with open(TEMP / "class_name", "r") as f:
        folder = f.read().replace("\n", "")
    return WORKFLOW + "resources/classes/" + folder + "/" + wildcards.file # DEBUG POSIX(WORKFLOW / "resources" / "classes" / folder / wildcards.file)


def class_files(wildcards):
    """

    """
    checkpoints.class_name.get(**wildcards)
    with open(TEMP / "class_name", "r") as f:
        folder = f.read().replace("\n", "")
    return [POSIX(TEX / file) for file in CLASSFILES.get(folder, [])]


checkpoint class_name:
    message:
        "Inferring document class..."
    output:
        temp(POSIX(TEMP / "class_name"))
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
        temp(POSIX(TEX / "{file}"))
    wildcard_constraints:
        file="|".join([file for files in CLASSFILES.values() for file in files])
    shell:
        "cp {input[0]} {output}"
