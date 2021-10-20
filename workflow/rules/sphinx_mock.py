"""
Mocks the imports of several variables / objects if the scripts
are imported during a Sphinx documentation build.

"""
try:

    assert __sphinx_docs_build__

    from pathlib import Path

    class workflow:
        class _showyourwork:
            snakefile = "/path/to/showyourwork/workflow"

        modules = {"showyourwork": _showyourwork}
        basedir = "/path/to/repo/"

    class abspaths:
        workflow = Path("showyourwork")

    class relpaths:
        tex = Path("src")
        temp = Path(".showyourwork")

    config = {}

    get_repo_url = lambda: "https://github.com/user/repo"

except:

    pass