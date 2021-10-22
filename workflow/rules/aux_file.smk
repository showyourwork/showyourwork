from pathlib import Path


localrules: aux_file


rule aux_file:
    """
    Copies auxiliary tex files to the ``src`` directory.
    
    """
    message:
        "Copying auxiliary tex file `{output}`..."
    input:
        lambda wildcards: posix(abspaths.workflow / "resources" / "tex" / Path(wildcards.file).name),
    wildcard_constraints:
        file="|".join(files.aux),
    output:
        temp("{file}"),
    shell:
        "cp {input} {output}"
