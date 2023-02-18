from showyourwork import test_util


def test_tex_build() -> None:
    test_util.run(
        "tests/projects/plugins/tex/build",
        snakemake_args=["ms.pdf"],
    )
