PLAN
----

Clone all user rules. Clones should be identical, except 

- Their names have a `syw__zenodo_` prefix.

- They precede the corresponding user rule in `ruleorder`.

- Their `script` points to a script that downloads the 
output from Zenodo (and `shell` and `notebook` are set to `None`).

- They have an extra `input` entry: a function that calls the
`dag` checkpoint (which blocks until the DAG is built) and returns
`[]` if and only if all of the rule's outputs exist on Zenodo and
the hash of the current rule matches the hash of the rule that generated
those outputs. Otherwise, returns a string corresponding to an 
unsatisfiable dependency (which causes the rule to be ignored).

- They have an extra `parameter` entry encoding all the information needed
to download the rule outputs from Zenodo in the `script`. This entry 
is a function that, like the input function above, blocks until the DAG 
is built, at which point we can calculate the rule hash and resolve the
outputs.

Then, to each user rule, add a dummy `input` entry that again blocks until
the DAG is built, and records the rule's resolved output and hash into global variables
if the job producing the output is actually run (we can query `dag.needrun(job)`
for this). Create a rule that is always run at the end of the workflow to
upload those outputs to Zenodo.



Also, figure out how to ensure the re-downloading of datasets if their Zenodo ID changes.