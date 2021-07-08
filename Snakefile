import glob
from pathlib import Path


TEXFILES = glob.glob(str(Path("tex") / "*"))
FIGFILES = glob.glob(str(Path("figures") / "*"))
TESTFILES = glob.glob(str(Path("tests") / "*"))
STYLEFILES = glob.glob(str(Path(".cortex") / "styles" / "*"))



rule meta:
    output: ".cortex/data/meta.json"
    shell: "cd .cortex && python user_prompt.py"

rule stylesheet:
    input: ancient(".cortex/data/meta.json")
    output: "tex/cortex.sty"
    shell: "cd .cortex && python generate_sty.py"

rule manuscript:
    input: ancient(".cortex/data/meta.json")
    output: "tex/ms.tex"
    shell: "cd .cortex && python generate_tex.py"

rule pdf:
    input: "tex/ms.tex", "tex/cortex.sty", TEXFILES, FIGFILES, TESTFILES
    output: "ms.pdf"
    run: 
        for f in STYLEFILES:
            shell("cp {} tex/".format(f))
        shell("tectonic -o . tex/ms.tex")
        for f in STYLEFILES:
            shell("rm tex/{}".format(Path(f).name))