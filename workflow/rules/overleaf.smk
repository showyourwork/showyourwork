from utils import paths


rule:
    """

    """
    name:
        "syw__overleaf_push"
    message:
        "Pushing file changes to Overleaf..."
    conda:
        "../envs/main.yml"
    script:
        "../scripts/overleaf_push.py"


rule:
    """

    """
    name:
        "syw__overleaf_pull"
    message:
        "Pulling file changes from Overleaf..."
    conda:
        "../envs/main.yml"
    script:
        "../scripts/overleaf_pull.py"