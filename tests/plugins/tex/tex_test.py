import pytest
from showyourwork.plugins.tex.tex import _add_to_preamble


@pytest.mark.parametrize("h", ["", "% a\n", "% a\n% b\n\n"])
@pytest.mark.parametrize("c", ["", "[a]", "[a,b]", "[a,b,c]"])
def test_add_to_preamble(h: str, c: str) -> None:
    value = _add_to_preamble(
        "PREAMBLE",
        r"""%s\documentclass%s{article}
\begin{document}
\end{document}
"""
        % (h, c),
    )
    expected = r"""%s\documentclass%s{article}

PREAMBLE

\begin{document}
\end{document}
""" % (
        h,
        c,
    )
    assert value == expected
