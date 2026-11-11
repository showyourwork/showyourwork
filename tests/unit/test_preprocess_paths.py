"""Copy commands for static figures must survive paths containing spaces."""

import importlib.util
import shlex
import subprocess

import pytest

from showyourwork import paths

SCRIPT = paths.showyourwork().workflow / "scripts" / "preprocess.py"

XML = """<FIGURE>
<GRAPHICS>figures/labeled.pdf</GRAPHICS>
<LABEL>fig:labeled</LABEL>
</FIGURE>
<GRAPHICS>figures/free.pdf</GRAPHICS>
"""


@pytest.fixture
def repo(tmp_path, monkeypatch):
    root = tmp_path / "my paper"
    for d in ("src/tex/figures", "src/static", ".showyourwork/preprocess"):
        (root / d).mkdir(parents=True)
    for name in ("labeled", "free"):
        (root / f"src/static/{name}.pdf").write_text(name)
    (root / ".showyourwork/preprocess/showyourwork.xml").write_text(XML)
    user = paths.user
    monkeypatch.setattr(paths, "user", lambda: user(root))
    return root


def test_static_copy_commands_run_in_repo_with_space(repo, monkeypatch):
    spec = importlib.util.spec_from_file_location("syw_preprocess", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(
        module,
        "config",
        {"dependencies": {}, "datasets": {}, "script_extensions": [], "scripts": {}},
        raising=False,
    )

    tree = module.get_json_tree(repo / ".showyourwork/preprocess/showyourwork.xml")

    commands = [
        tree["figures"]["fig:labeled"]["command"],
        tree["figures"]["free-floating-static"]["command"],
    ]
    for command in commands:
        assert shlex.split(command)[0] == "cp"
        subprocess.run(command, shell=True, cwd=repo, check=True)
    assert (repo / "src/tex/figures/labeled.pdf").read_text() == "labeled"
    assert (repo / "src/tex/figures/free.pdf").read_text() == "free"
