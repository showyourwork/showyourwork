import subprocess
from pathlib import Path
import os

# Top level
ROOT = Path(__file__).parents[1].absolute()

# Generate the DAG. This is in general a terrible idea -- we shouldn't
# be calling Snakemake within a Snakemake build! But hey, it works.
dag = subprocess.check_output(["snakemake", "--dag"], cwd=ROOT).decode()

# Trim the file names to reduce the image width
dag = dag.replace("figure: figures/", "")
dag = dag.replace("test: tests/", "")

# Save to a temporary file
with open("dag.dag", "w") as f:
    print(dag, file=f)

# Generate the PDF
with open("dag.pdf", "wb") as f:
    f.write(subprocess.check_output(["dot", "-Tpdf", "dag.dag"]))

# Remove the tempfile
os.remove("dag.dag")
