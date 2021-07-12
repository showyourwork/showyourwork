import glob
from pathlib import Path
import os


# Debug mode?
debug = bool(int(os.environ.get("CTX_DEBUG", "0")))

# Paths
TEX_FILES = glob.glob(str(Path("tex") / "*"))
FIGURE_SCRIPTS = glob.glob(str(Path("figures") / "*.py"))
PYTHON_TESTS = glob.glob(str(Path("tests") / "*"))


# Generate the metadata file from direct user input
rule meta:
    output: ".cortex/data/meta.json"
    shell: "cd .cortex && python user_prompt.py"

# Generate figures (TODO)
rule figures:
    input: "figures/{figure}.py"
    output: "figures/{figure}{suffix,(_.*)?}.pdf"
    shell: "cd figures && python {wildcards.figure}.py"

# Compile the PDF
rule pdf:
    input: "tex/ms.tex", "tex/bib.bib", TEX_FILES, FIGURE_SCRIPTS, PYTHON_TESTS
    output: "ms.pdf"
    run: 
        shell("cd .cortex && python generate_sty.py")
        shell("cd .cortex && python check_git_status.py")
        STYLE_FILES = glob.glob(str(Path(".cortex") / "styles" / "*"))
        TEMP_FILES = glob.glob(str(Path("figures") / "*.py.cortex"))

        def clean():
            if not debug:
                for f in STYLE_FILES:
                    shell("rm tex/{}".format(Path(f).name))
                for f in TEMP_FILES:
                    shell("rm {}".format(f))

        try:
            for f in STYLE_FILES:
                shell("cp {} tex/".format(f))
            if debug:
                shell("tectonic --print -o . tex/ms.tex")
            else:
                shell("tectonic -o . tex/ms.tex")
            clean()
        except Exception as e:
            clean()
            raise e