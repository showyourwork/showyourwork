"""
Defines the rule ``syw__render_dag`` to render the article DAG.

Runs the script :doc:`render_dag` to generate ``dag.pdf``.

"""
from showyourwork import paths

# All figure output files
figures = config["tree"]["figures"].keys()

rule:
    """
    Render the article DAG.

    """
    name:
        "syw__render_dag"
    message:
        "Rendering the article DAG..."
    input:
        config["ms_tex"],
        config["dependencies"][config["ms_tex"]],
        WORKFLOW_GRAPH,
        "showyourwork.yml"
    output:
       "dag.pdf"
    script:
        "../scripts/render_dag.py"
