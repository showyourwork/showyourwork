"""
Sets up the internal conda environment.

"""
from showyourwork import paths


rule:
    """
    Sets up the internal conda environment.

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
