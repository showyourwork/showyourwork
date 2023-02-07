from pathlib import Path

from showyourwork import paths

repo_path = paths.repo(config).root
dynamic_path = paths.work(config).root
build_path = paths.work(config).build

# Default script rules
scripts = {
    "py": "python {script} {output}",
    "ipynb": "jupyter execute {script}",
}
scripts = dict(scripts, **config.get("scripts", {}))

for dynamic in config.get("dynamic", []):
    script = repo_path / dynamic["script"]
    name = f"syw__dynamic_{paths.path_to_rule_name(dynamic['script'])}"
    input = [repo_path / f for f in dynamic.get("input", [])],
    output = dynamic.get("output", [])
    if isinstance(output, str):
        output = [repo_path / dynamic["output"]]
    else:
        output = [repo_path / f for f in dynamic.get("output", [])]
    conda = repo_path / dynamic.get("conda", config.get("conda", None))

    suffix = Path(script).suffix[1:]
    command = scripts[suffix].format(script=script, output=" ".join(map(str, output)))

    if conda is None:
        rule:
            name:
                name
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
            input:
                input, script=script
            output:
                output
            conda:
                conda
            shell:
                command
