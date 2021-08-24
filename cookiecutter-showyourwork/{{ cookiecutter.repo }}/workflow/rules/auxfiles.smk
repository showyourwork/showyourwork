# Auxiliary files we copy over to the user tex/ directory
aux_files = [posix(TEX / file.name) for file in (WORKFLOW / "resources" / "tex").glob("*.*")]


rule aux_file:
    message:
        "Copying auxiliary tex file `{output}`..."
    input:
        lambda wildcards: posix(WORKFLOW / "resources" / wildcards.file)
    wildcard_constraints:
        file="|".join(aux_files)
    output:
        temp("{file}")
    shell:
        "cp {input} {output}"
