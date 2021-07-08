import glob
from pathlib import Path


TEXFILES = glob.glob(str(Path("tex") / "*"))
FIGFILES = glob.glob(str(Path("figures") / "*"))
TESTFILES = glob.glob(str(Path("tests") / "*"))
STYLEFILES = glob.glob(str(Path(".cortex") / "styles" / "*"))

rule meta:
    output: "data/meta.json"
    shell: "python user_prompt.py"

rule stylesheet:
    input: ancient("data/meta.json")
    output: "tex/cortex.sty"
    shell: "python generate_sty.py"

rule manuscript:
    input: ancient("data/meta.json")
    output: "tex/ms.tex", "tex/cortex.sty"
    shell: "python generate_tex.py"

rule pdf:
    input: TEXFILES, FIGFILES, TESTFILES
    output: "ms.pdf"
    run: 
        for f in STYLEFILES:
            shell("cp {} tex/".format(f))
        shell("tectonic -o . tex/ms.tex")
        for f in STYLEFILES:
            shell("rm tex/{}".format(Path(f).name))