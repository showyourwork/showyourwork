import glob
from pathlib import Path


TEXFILES = glob.glob(str(Path("tex") / "*"))
FIGFILES = glob.glob(str(Path("figures") / "*"))
TESTFILES = glob.glob(str(Path("tests") / "*"))
STYLEFILES = glob.glob(str(Path(".cortex") / "styles" / "*"))



# Generate the metadata file from direct user input
rule meta:
    output: ".cortex/data/meta.json"
    shell: "cd .cortex && python user_prompt.py"

# Generate the `cortex.sty` stylesheet
rule stylesheet:
    input: ancient(".cortex/data/meta.json")
    output: "tex/cortex.sty"
    shell: "cd .cortex && python generate_sty.py"

# Generate the skeleton `ms.tex` manuscript
rule manuscript:
    input: ancient(".cortex/data/meta.json")
    output: "tex/ms.tex"
    shell: "cd .cortex && python generate_tex.py"

# Compile the PDF
rule pdf:
    input: "tex/ms.tex", "tex/cortex.sty", TEXFILES, FIGFILES, TESTFILES
    output: "ms.pdf"
    run: 
        for f in STYLEFILES:
            shell("cp {} tex/".format(f))
        shell("tectonic -o . tex/ms.tex")
        for f in STYLEFILES:
            shell("rm tex/{}".format(Path(f).name))