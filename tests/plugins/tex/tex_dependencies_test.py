import json
from contextlib import contextmanager
from pathlib import Path
from typing import Generator

from showyourwork import paths, test_util

plugin_id = "showyourwork.plugins.tex"


@contextmanager
def generate_dependencies(contents: str) -> Generator[Path, None, None]:
    with test_util.temporary_project() as d:
        p = paths.work({"working_directory": d})
        open(p.manuscript, "w").write(contents)
        test_util.run_snakemake(
            str(paths.package_data(plugin_id, "workflow", "Snakefile")),
            ["_syw_plug_tex_dependencies", "--config", f"working_directory={d}"],
            cwd=d,
        )
        assert (p.plugin(plugin_id, "xml") / "showyourwork.xml").is_file()
        json_file = p.dependencies
        assert json_file.is_file()
        yield json_file


def test_basic() -> None:
    with generate_dependencies(
        r"""
\documentclass{article}
\usepackage{showyourwork}
\begin{document}
Test.
\end{document}
"""
    ):
        pass


def test_unlabeled_figure() -> None:
    with generate_dependencies(
        r"""
\documentclass{article}
\usepackage{showyourwork}
\begin{document}
Test.
\begin{figure}
\includegraphics{test.png}
\end{figure}
\end{document}
            """
    ) as json_file:
        with open(json_file, "r") as f:
            data = json.load(f)
        assert data["figures"] == {}
        assert data["files"] == []
        assert len(data["unlabled"]) == 1
        assert Path(data["unlabled"][0]).name == "test.png"


def test_labeled_figure() -> None:
    with generate_dependencies(
        r"""
\documentclass{article}
\usepackage{showyourwork}
\begin{document}
Test.
\begin{figure}
\includegraphics{test.png}
\label{fig:test}
\end{figure}
\end{document}
            """
    ) as json_file:
        with open(json_file, "r") as f:
            data = json.load(f)
        assert "fig:test" in data["figures"]
        assert len(data["figures"]["fig:test"]) == 1
        assert Path(data["figures"]["fig:test"][0]).name == "test.png"
        assert data["files"] == []
        assert data["unlabled"] == []
