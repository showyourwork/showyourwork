repo_path = SYW__REPO_PATHS.root
dynamic_path = SYW__WORK_PATHS.root
build_path = SYW__WORK_PATHS.build

# Default script rules
scripts = {
    "py": "python {script} {output}",
    "ipynb": "jupyter execute {script}",
}
scripts = dict(scripts, **config.get("scripts", {}))

for dynamic in config.get("dynamic", []):
    script = repo_path / dynamic["script"]
    name = utils.rule_name("dynamic", document=dynamic["script"])
    message = f"Running script '{dynamic['script']}'"
    input = [repo_path / f for f in dynamic.get("input", [])],
    output = dynamic.get("output", [])
    if isinstance(output, str):
        output = [output]
    conda = repo_path / dynamic.get("conda", config.get("conda", None))

    suffix = Path(script).suffix[1:]
    command = scripts[suffix].format(script=script, output=" ".join(map(str, output)))

    if conda is None:
        rule:
            name:
                name
            message:
                message
            input:
                input, script=script
            output:
                output
            shell:
                command

    else:
        rule:
            name:
                name
            message:
                message
            input:
                input, script=script
            output:
                output
            conda:
                conda
            shell:
                command
