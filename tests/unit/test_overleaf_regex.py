import re

import pytest

from showyourwork.overleaf import (
    OVERLEAF_BLANK_PROJECT,
    OVERLEAF_BLANK_PROJECT_REGEX_TEMPLATE,
    clone,
    get_remote_branch,
)


def test_current_template():
    assert re.match(OVERLEAF_BLANK_PROJECT_REGEX_TEMPLATE, OVERLEAF_BLANK_PROJECT)


def test_old_template():
    assert re.match(
        OVERLEAF_BLANK_PROJECT_REGEX_TEMPLATE,
        r"""\documentclass{article}
\usepackage{graphicx} % Required for inserting images

\title{showyourwork test}
\author{Dan Foreman-Mackey}
\date{March 2023}

\begin{document}

\maketitle

\section{Introduction}

\end{document}
""",
    )


@pytest.mark.parametrize(
    "stdout, expected",
    [
        ("ref: refs/heads/main HEAD\nabc123\tHEAD\n", "main"),
        ("ref: refs/heads/master HEAD\nabc123\tHEAD\n", "master"),
        ("abc123\tHEAD\n", "master"),
    ],
)
def test_get_remote_branch(monkeypatch, stdout, expected):
    def fake_get_stdout(args, cwd=None, secrets=(), callback=None, **kwargs):
        assert args == [
            "git",
            "ls-remote",
            "--symref",
            "https://git:token@git.overleaf.com/project",
            "HEAD",
        ]
        return callback(0, stdout, "")

    monkeypatch.setattr("showyourwork.overleaf.get_stdout", fake_get_stdout)

    assert (
        get_remote_branch("https://git:token@git.overleaf.com/project", "token")
        == expected
    )


def test_clone_uses_remote_default_branch(monkeypatch, tmp_path):
    calls = []

    def fake_get_stdout(args, cwd=None, secrets=(), callback=None, **kwargs):
        calls.append(args)

        if args[:3] == ["git", "ls-remote", "--symref"]:
            return callback(0, "ref: refs/heads/main HEAD\nabc123\tHEAD\n", "")
        elif args == ["git", "init"]:
            return ""
        elif args == ["git", "checkout", "-b", "main"]:
            return ""
        elif args == ["git", "pull", "https://git:token@git.overleaf.com/project"]:
            return callback(0, "", "")
        else:
            raise AssertionError(args)

    monkeypatch.setattr("showyourwork.overleaf.get_stdout", fake_get_stdout)
    monkeypatch.setenv("OVERLEAF_TOKEN", "token")

    branch = clone("project", path=tmp_path)

    assert branch == "main"
    assert ["git", "checkout", "-b", "main"] in calls
    assert ["git", "pull", "https://git:token@git.overleaf.com/project"] in calls
