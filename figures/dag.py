import subprocess
from pathlib import Path

ROOT = Path(__file__).parents[1].absolute()
ps = subprocess.Popen(("snakemake", "--dag"), stdout=subprocess.PIPE, cwd=ROOT)
with open("dag.pdf", "wb") as f:
    f.write(subprocess.check_output(("dot", "-Tpdf"), stdin=ps.stdout))
ps.wait()
