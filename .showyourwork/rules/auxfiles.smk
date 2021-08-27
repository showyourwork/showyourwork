from pathlib import Path


localrules:
    aux_file,


rule aux_file:
    message:
        "Copying auxiliary tex file `{output}`..."
    input:
        lambda wildcards: POSIX(WORKFLOW / "resources" / "tex" / Path(wildcards.file).name),
    wildcard_constraints:
        file="|".join(AUXFILES),
    output:
        temp("{file}"),
    shell:
        "cp {input} {output}"
