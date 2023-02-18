from showyourwork import test_util


def test_render() -> None:
    test_util.run("tests/projects/render", snakemake_args=["syw__render_dag"])
