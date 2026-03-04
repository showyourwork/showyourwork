"""Integration tests for figures with different placement specifiers.

Regression tests for issue #626: the float package's [H] specifier
internally replaces \\endfigure, which bypassed the XML logging in
preprocess.tex. The fix uses etoolbox environment hooks instead.
"""

from helpers import (
    ShowyourworkRepositoryActions,
    TemporaryShowyourworkRepository,
)


class TestFigurePlacementH(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """Test that figures with [H] placement build without XML parsing errors."""

    # No need to test this on CI
    local_build_only = True

    def customize(self):
        # Add float package to preamble
        ms = self.cwd / "src" / "tex" / "ms.tex"
        content = ms.read_text()
        content = content.replace(
            r"\usepackage{showyourwork}",
            r"\usepackage{showyourwork}" + "\n" + r"\usepackage{float}",
        )
        ms.write_text(content)

        # Add a figure with [H] placement
        lines = ms.read_text().splitlines()
        figure_block = "\n".join(
            [
                r"\begin{figure}[H]",
                r"\script{test_figure.py}",
                r"\begin{centering}",
                r"\includegraphics[width=\linewidth]{figures/test_figure.pdf}",
                r"\caption{A test figure with [H] placement.}",
                r"\label{fig:h_placement}",
                r"\end{centering}",
                r"\end{figure}",
                "",
            ]
        )

        for idx, line in enumerate(lines):
            if line.strip() == r"\end{document}":
                lines.insert(idx, figure_block)
                break

        ms.write_text("\n".join(lines) + "\n")

        # Add the figure script
        self.add_figure_script()


class TestFigurePlacementMixed(
    TemporaryShowyourworkRepository, ShowyourworkRepositoryActions
):
    """Test figures with mixed placement specifiers including [H]."""

    # No need to test this on CI
    local_build_only = True

    def customize(self):
        # Add float package to preamble
        ms = self.cwd / "src" / "tex" / "ms.tex"
        content = ms.read_text()
        content = content.replace(
            r"\usepackage{showyourwork}",
            r"\usepackage{showyourwork}" + "\n" + r"\usepackage{float}",
        )
        ms.write_text(content)

        # Add figures with different placement specifiers
        lines = ms.read_text().splitlines()

        # Figure without placement specifier (default)
        fig1 = "\n".join(
            [
                r"\begin{figure}",
                r"\script{test_figure.py}",
                r"\begin{centering}",
                r"\includegraphics[width=\linewidth]{figures/test_figure.pdf}",
                r"\caption{Figure with default placement.}",
                r"\label{fig:default}",
                r"\end{centering}",
                r"\end{figure}",
                "",
            ]
        )

        # Figure with [ht!]
        fig2 = "\n".join(
            [
                r"\begin{figure}[ht!]",
                r"\script{test_figure.py}",
                r"\begin{centering}",
                r"\includegraphics[width=\linewidth]{figures/test_figure.pdf}",
                r"\caption{Figure with [ht!] placement.}",
                r"\label{fig:ht}",
                r"\end{centering}",
                r"\end{figure}",
                "",
            ]
        )

        # Figure with [H] (float package)
        fig3 = "\n".join(
            [
                r"\begin{figure}[H]",
                r"\script{test_figure.py}",
                r"\begin{centering}",
                r"\includegraphics[width=\linewidth]{figures/test_figure.pdf}",
                r"\caption{Figure with [H] placement.}",
                r"\label{fig:h_exact}",
                r"\end{centering}",
                r"\end{figure}",
                "",
            ]
        )

        for idx, line in enumerate(lines):
            if line.strip() == r"\end{document}":
                lines.insert(idx, fig3)
                lines.insert(idx, fig2)
                lines.insert(idx, fig1)
                break

        ms.write_text("\n".join(lines) + "\n")

        # Add the figure script
        self.add_figure_script()
