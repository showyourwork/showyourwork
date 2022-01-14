from pathlib import Path


workflow = Path(__file__).absolute().parents[1]
rules = workflow / "rules"
checkpoints = workflow / "checkpoints"
resources = workflow / "resources"

showyourwork = workflow.parents[0]

user = showyourwork.parents[0]
src = user / "src"
tex = src / "tex"
figure_scripts = src / "figure-scripts"
static_figures = src / "static-figures"
figures = tex / "figures"


temp = user / ".showyourwork"
temp.mkdir(exist_ok=True)

preprocess = temp / "preprocess"
preprocess.mkdir(exist_ok=True)

compile = temp / "compile"
compile.mkdir(exist_ok=True)