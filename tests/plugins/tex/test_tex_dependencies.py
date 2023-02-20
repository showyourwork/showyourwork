from showyourwork import test_util


def test_tex_dependencies() -> None:
    test_util.run(
        "tests/projects/plugins/tex/dependencies",
        snakemake_args=["sywplug_tex__dependencies"],
    )
