
from pathlib import Path
import subprocess

envfile = snakemake.params.envfile
envdir = snakemake.params.envdir
output = snakemake.output[0]

with open(envfile, "r") as f:
    envfile = f.read()

for file in Path(envdir).glob("*.yaml"):
    with open(file, "r") as f:
        envfile_ = f.read()
    if envfile == envfile_:
        env = str(file.parents[0] / file.stem)
        conda_prefix = subprocess.run(["conda", "info", "--base"], stdout=subprocess.PIPE).stdout.decode().replace("\n", "")
        conda_activate = f". {conda_prefix}/etc/profile.d/conda.sh && conda activate {env} && "
        break
else:
    conda_activate = ""

with open(output, "w") as f:
    f.write(conda_activate)