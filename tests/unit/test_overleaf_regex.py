import re

from showyourwork.overleaf import (
    OVERLEAF_BLANK_PROJECT,
    OVERLEAF_BLANK_PROJECT_REGEX_TEMPLATE,
)


def test_current_template():
    assert re.match(
        OVERLEAF_BLANK_PROJECT_REGEX_TEMPLATE, OVERLEAF_BLANK_PROJECT
    )


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
