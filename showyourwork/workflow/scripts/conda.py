"""
Figure out which ``conda`` environment in the ``.snakemake`` folder corresponds
to ``showyourwork/showyourwork/workflow/envs/environment.yml`` and record the
command to activate it in a temporary file. This allows us to temporarily
activate it to run commands like ``tectonic``, ``dot`` (from ``graphviz``), and
``convert`` (from ``imagemagick``).

.. todo::

    This is hacky, since we should just be running these commands
    within isolated rules that explicitly depend on that environment! This is
    tricky because currently we're running things like ``tectonic`` deep in
    the middle of certain rules that also require the ``showyourwork`` module.
    We'll have to divide those rules (like ``compile``) into several sub-rules
    so we can compartmentalize the build and take advantage of this custom
    ``conda`` environment.

"""

import subprocess
from pathlib import Path

if __name__ == "__main__":

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
            conda_prefix = (
                subprocess.run(
                    ["conda", "info", "--base"], stdout=subprocess.PIPE
                )
                .stdout.decode()
                .replace("\n", "")
            )
            conda_activate = f'. {conda_prefix}/etc/profile.d/conda.sh && conda activate "{env}" && '
            break
    else:
        conda_activate = ""

    with open(output, "w") as f:
        f.write(conda_activate)
