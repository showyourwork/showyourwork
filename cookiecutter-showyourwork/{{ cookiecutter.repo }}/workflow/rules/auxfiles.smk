rule aux_file:
    message:
        "Copying auxiliary tex file `{output}`..."
    input:
        lambda wildcards: WORKFLOW + "resources/" + wildcards.file # DEBUG POSIX(WORKFLOW / "resources" / wildcards.file)
    wildcard_constraints:
        file="|".join(AUXFILES)
    output:
        temp("{file}")
    shell:
        "cp {input} {output}"
