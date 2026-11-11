"""Tests that edits to the manuscript are picked up by a rebuild."""

from helpers import TemporaryShowyourworkRepository


class TestEditManuscript(TemporaryShowyourworkRepository):
    """Editing ``ms.tex`` after a build must recompile the article (#778)."""

    local_build_only = True

    def check_build(self):
        ms_tex = self.cwd / "src" / "tex" / "ms.tex"
        compiled_tex = self.cwd / ".showyourwork" / "compile" / "ms.tex"
        ms_pdf = self.cwd / "ms.pdf"
        assert ms_pdf.exists()
        pdf_mtime = ms_pdf.stat().st_mtime

        ms_tex.write_text(ms_tex.read_text() + "\nSentinel text xyzzy.\n")
        self.git_commit()
        self.build_local()

        assert "xyzzy" in compiled_tex.read_text()
        assert ms_pdf.stat().st_mtime > pdf_mtime
