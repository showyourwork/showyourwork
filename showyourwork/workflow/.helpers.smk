import showyourwork
from showyourwork.utils import glob
from showyourwork.constants import TEMP, USER, FIGURE_EXTENSIONS
from pathlib import Path

# User config
verbose = (config.get("verbose", "false").lower() == "true")


def figure_scripts():
    return showyourwork.meta.get_script_metadata(clobber=False, verbose=verbose)[
        "figures"
    ].values()


def figures(wildcards):
    checkpoints.script_info.get(**wildcards)
    figures = []
    for entry in figure_scripts():
        figures += entry["files"]
    return figures


def figure_script(wildcards):
    checkpoints.script_info.get(**wildcards)
    figure = wildcards.figure
    for entry in figure_scripts():
        if figure in entry["files"]:
            return entry["script"]
    raise ValueError(
        "Input script not found for output figure `{}`.".format(figure)
    )


def figure_data(wildcards):
    script = figure_script(wildcards)
    for prop_name in dir(rules):
        if not prop_name.startswith("__"):
            prop = getattr(rules, prop_name) #.rule
            if hasattr(prop, "rule"):
                syw_rule = prop.rule
                if Path(syw_rule.snakefile).name == "Snakefile":
                    data_files = syw_rule.output
                    if script in syw_rule.params.get("scripts", []):
                        return data_files
    return []


def cache_cmd(wildcards, input, output):
    other = []
    for entry in figure_scripts():
        if entry["script"] == input[0]:
            other = [
                Path(file).name for file in entry["files"] if file != output[0]
            ]
    if len(other):
        return " && ".join(
            [f"mv {file} {TEMP / 'figures'}" for file in other]
        )
    else:
        return ":"


def cached_figure(wildcards, input, output):
    return str(TEMP / output[0])


def script_name(wildcards, input):
    return Path(input[0]).name


def run_pdf():
    showyourwork.utils.make_pdf(
        tmpdir=TEMP / "tex",
        verbose=verbose,
        publish=True,
        **showyourwork.meta.get_metadata(clobber=False),
    )


def run_repo_info():
    showyourwork.meta.get_repo_metadata(clobber=True)


def run_script_info():
    showyourwork.meta.get_script_metadata(clobber=True, verbose=verbose)


def run_metadata():
    showyourwork.meta.get_metadata(clobber=True)
