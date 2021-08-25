import jinja2
import json


# Params defined in `../rules/pdf.smk`
WORKFLOW = snakemake.params["WORKFLOW"]
TEMP = snakemake.params["TEMP"]
TEX = snakemake.params["TEX"]


# Generate the stylesheet
env = jinja2.Environment(
    block_start_string="((*",
    block_end_string="*))",
    variable_start_string="((-",
    variable_end_string="-))",
    comment_start_string="((=",
    comment_end_string="=))",
    trim_blocks=True,
    autoescape=False,
    loader=jinja2.FileSystemLoader(WORKFLOW / "resources" / "templates"),
)
with open(TEMP / "meta.json", "r") as f:
    jinja_kwargs = json.load(f)
with open(TEX / "showyourwork.sty", "w") as f:
    sty = env.get_template("showyourwork.sty").render(**jinja_kwargs)
    print(sty, file=f)
