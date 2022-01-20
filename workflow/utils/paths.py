from pathlib import Path
from os.path import realpath


# Workflow paths
workflow = Path(realpath(__file__)).absolute().parents[1]
rules = workflow / "rules"
checkpoints = workflow / "checkpoints"
resources = workflow / "resources"


# The showyourwork submodule
showyourwork = workflow.parents[0]


# User paths
user = showyourwork.parents[0]
src = user / "src"
tex = src / "tex"
data = src / "data"
figure_scripts = src / "figure-scripts"
static_figures = src / "static-figures"
figures = tex / "figures"


# Temporary paths
temp = user / ".showyourwork"
temp.mkdir(exist_ok=True)
preprocess = temp / "preprocess"
preprocess.mkdir(exist_ok=True)
compile = temp / "compile"
compile.mkdir(exist_ok=True)
logs = temp / "logs"
logs.mkdir(exist_ok=True)
zenodo = temp / "zenodo"
zenodo.mkdir(exist_ok=True)
zenodo_sandbox = temp / "zenodo_sandbox"
zenodo_sandbox.mkdir(exist_ok=True)