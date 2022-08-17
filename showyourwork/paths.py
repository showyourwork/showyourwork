from os.path import realpath
from pathlib import Path

from . import git


# showyourwork paths
class showyourwork:
    """
    Paths to various directories within the ``showyourwork`` package.

    """

    def __init__(self):
        self.module = Path(realpath(__file__)).absolute().parents[0]
        self.workflow = self.module / "workflow"
        self.rules = self.workflow / "rules"
        self.resources = self.workflow / "resources"
        self.envs = self.workflow / "envs"
        self.cookiecutter = self.module / "cookiecutter-showyourwork"


# User paths
class user:
    """
    Paths to various directories within the user's repository.

    """

    def __init__(self, path=None):
        """

        Args:
            path (str): The path to the top level of the user's repository
                (if running outside of the repository).
        """
        if path is None:
            root = git.get_repo_root()
            if root == "unknown":
                raise Exception("Not in a git repo.")
            else:
                path = Path(root).absolute()

        # Repo paths
        self.repo = Path(path)
        self.src = self.repo / "src"
        self.tex = self.src / "tex"
        self.data = self.src / "data"
        self.scripts = self.src / "scripts"
        self.static = self.src / "static"
        self.figures = self.tex / "figures"
        self.output = self.tex / "output"

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
        self.sandbox = self.temp / "sandbox"
        self.sandbox.mkdir(exist_ok=True)
        self.sandbox_ids = self.sandbox / "ids"
        self.sandbox_ids.mkdir(exist_ok=True)
        self.overleaf = self.temp / "overleaf"
        self.overleaf.mkdir(exist_ok=True)
        self.flags = self.temp / "flags"
        self.flags.mkdir(exist_ok=True)
