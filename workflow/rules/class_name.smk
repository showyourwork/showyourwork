import os
import re


localrules: class_name


checkpoint class_name:
    """
    Infers the name of the LaTeX document class for the article
    by doing a regex on the TeX file. This helps ``showyourwork``
    decide which class files to include in the build.

    """
    message:
        "Inferring document class..."
    output:
        posix(relpaths.temp / "class_name"),
    priority: 100
    run:
        with open(relpaths.tex / config["ms_name"], "r") as f:
            lines = f.readlines()
            for line in lines:
                match = re.match(
                    "[ \t]*\\\documentclass\[?.*?\]?\{(.*?)\}", line
                )
                if hasattr(match, "groups"):
                    name = match.groups()[0]
                    break
            else:
                raise ValueError(
                    "Unable to determine document class in `{}`.".format(config["ms_name"])
                )
            if not relpaths.temp.exists():
                os.mkdir(relpaths.temp)
            with open(relpaths.temp / "class_name", "w") as f:
                print(name, file=f)



