import pytest

from showyourwork.cli.conda_env import parse_syw_spec


@pytest.mark.parametrize(
    "spec,expected",
    [
        (None, "showyourwork"),
        ("0.3.0", "showyourwork==0.3.0"),
        (
            "893eda2",
            "git+https://github.com/showyourwork/showyourwork.git@893eda2#egg=showyourwork",
        ),
        (
            "893eda204de56455002e5ef2d1d49e0b7349d37a",
            "git+https://github.com/showyourwork/showyourwork.git@893eda204de56455002e5ef2d1d49e0b7349d37a#egg=showyourwork",
        ),
        (
            "main",
            "git+https://github.com/showyourwork/showyourwork.git@main#egg=showyourwork",
        ),
        ({"pip": "0.3.0"}, "showyourwork==0.3.0"),
        (
            {"branch": "branch"},
            "git+https://github.com/showyourwork/showyourwork.git@branch#egg=showyourwork",
        ),
        (
            {"fork": "https://github.com/dfm/showyourwork"},
            "git+https://github.com/dfm/showyourwork#egg=showyourwork",
        ),
        (
            {"sha": "893eda2"},
            "git+https://github.com/showyourwork/showyourwork.git@893eda2#egg=showyourwork",
        ),
        (
            {"tag": "v0.3.0"},
            "git+https://github.com/showyourwork/showyourwork.git@v0.3.0#egg=showyourwork",
        ),
        (
            {
                "fork": "https://github.com/dfm/showyourwork",
                "branch": "branch",
            },
            "git+https://github.com/dfm/showyourwork@branch#egg=showyourwork",
        ),
    ],
)
def test_parse_syw_spec(spec, expected):
    assert parse_syw_spec(spec) == expected
