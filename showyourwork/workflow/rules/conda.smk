"""
Sets up the ``conda`` environment used to run certain showyourwork rules that
depend on non-python packages such as ``tectonic`` and ``graphviz``.

"""
from showyourwork import paths


rule:
    """
    Set up the ``conda`` environment used to run certain showyourwork rules that
    depend on non-python packages such as ``tectonic`` and ``graphviz``.
    This rule forces ``Snakemake`` to create the necessary ``conda`` environment
    and writes the command to activate it, e.g.,

    .. code-block::bash

        . <conda-path>/etc/profile.d/conda.sh && conda activate <env-path>

    in the temporary file ``.showyourwork/flags/SYW_CONDA``. We can then
    easily run commands like ``tectonic`` and ``dot`` within that environment
    within rules that depend on ``.showyourwork/flags/SYW_CONDA`` as an input.

    .. todo::

        This is hacky, since we should just be running these commands
        within isolated rules that explicitly depend on that environment! This is
        tricky because currently we're running things like ``tectonic`` deep in
        the middle of certain rules that also require the ``showyourwork`` module.
        We'll have to divide those rules (like ``compile``) into several sub-rules
        so we can compartmentalize the build and take advantage of this custom
        ``conda`` environment.

    """
    name:
        "syw__conda"
    message:
        "Setting up the internal conda environment..."
    output:
        paths.user().flags / "SYW__CONDA"
    conda:
        (paths.showyourwork().envs / "environment.yml").as_posix()
    params:
        envfile=(paths.showyourwork().envs / "environment.yml").as_posix(),
        envdir=(paths.user().snakemake / "conda").as_posix(),
    script:
        "../scripts/conda.py"
