localrules: aux_file

rule aux_file:
    message:
        "Copying auxiliary tex file `{output}`..."
    input:
        lambda wildcards: POSIX(WORKFLOW / "resources" / wildcards.file)
    wildcard_constraints:
        file="|".join(AUXFILES)
    output:
        temp("{file}")
    shell:
        "cp {input} {output}"
