from . import git
from pathlib import Path
from os.path import realpath


# showyourwork paths
class showyourwork:
    def __init__(self):
        self.module = Path(realpath(__file__)).absolute().parents[0]
        self.workflow = self.module / "workflow"
        self.rules = self.workflow / "rules"
        self.checkpoints = self.workflow / "checkpoints"
        self.resources = self.workflow / "resources"
        self.envs = self.workflow / "envs"
        self.cookiecutter = self.module / "cookiecutter-showyourwork"


# User paths
class user:
    def __init__(self, path=None):

        if path is None:
            root = git.get_repo_root()
            if root == "None":
                raise Exception("Not in a git repo.")
            else:
                path = Path(root).absolute()

        # Repo paths
        self.repo = Path(path)
        self.src = self.repo / "src"
        self.tex = self.src / "tex"
        self.data = self.src / "data"
        self.figure_scripts = self.src / "figures"
        self.static_figures = self.src / "static"
        self.figures = self.tex / "figures"

        # User home temp (for all repos)
        self.home_temp = Path.home() / ".showyourwork"
        self.home_temp.mkdir(exist_ok=True)
        self.env = self.home_temp / "env"

        # Temporary paths
        self.temp = self.repo / ".showyourwork"
        self.temp.mkdir(exist_ok=True)
        self.cache = self.temp / "cache"
        self.cache.mkdir(exist_ok=True)
        self.preprocess = self.temp / "preprocess"
        self.preprocess.mkdir(exist_ok=True)
        self.compile = self.temp / "compile"
        self.compile.mkdir(exist_ok=True)
        self.logs = self.temp / "logs"
        self.logs.mkdir(exist_ok=True)
        self.zenodo = self.temp / "zenodo"
        self.zenodo.mkdir(exist_ok=True)
        self.zenodo_ids = self.zenodo / "ids"
        self.zenodo_ids.mkdir(exist_ok=True)
        self.zenodo_sandbox = self.temp / "zenodo_sandbox"
        self.zenodo_sandbox.mkdir(exist_ok=True)
        self.zenodo_sandbox_ids = self.zenodo_sandbox / "ids"
        self.zenodo_sandbox_ids.mkdir(exist_ok=True)
        self.overleaf = self.temp / "overleaf"
        self.overleaf.mkdir(exist_ok=True)
        self.checkpoints = self.temp / "checkpoints"
        self.checkpoints.mkdir(exist_ok=True)