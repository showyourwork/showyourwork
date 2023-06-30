from showyourwork.testing import run_showyourwork


def test_tex_dependencies() -> None:
    run_showyourwork(
        "tests/projects/plugins/tex/dependencies",
        "sywplug_tex__dependencies",
        show_diff=True,
    )
