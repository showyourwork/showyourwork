"""
Defines the ``syw__fallback__`` rules.

These are fallback rules for manuscript dependencies, which raise an error
if they ever get invoked. Implementing these is necessary because of certain
edge cases involving the behavior of Snakemake when dependencies aren't
present.

For example, suppose the user runs ``showyourwork`` and successfully builds
their article, ``ms.pdf``. They then add a figure environment to the ``ms.tex``
manuscript with calls to ``\includegraphics{figure.pdf}`` and
``\script{figure.py}``. This tells ``showyourwork`` that (1) ``figure.pdf`` is
a build dependency of the article and (2) ``figure.py`` can be executed to
produce ``figure.pdf``. Now suppose the user forgets to create the file
``figure.py`` and executes ``showyourwork`` again. The workflow *should* fail,
but instead it completes successfully because Snakemake assures us there
is "Nothing to be done".

Why is that? Snakemake is aware that ``ms.tex`` is modified, so it attempts to
run the rule ``syw__compile`` to produce ``ms.pdf`` (the downstream
output of the modified ``ms.tex`` file). But that rule is ineligible for the
workflow run because one of its other dependencies (``figure.py``) does not
exist. The default behavior of Snakemake in this case is to skip this rule 
without throwing an error, since there may be other rules that can still
complete the job successfully. But it finds none, so the workflow simply completes
with the message "Nothing to be done".

This may be the intended behavior in some cases (I'm not sure!), but in our
case it's not: the workflow should fail and inform the user that there is no
rule to generate a necessary dependency of the article PDF.

So, long story short, this Snakefile implements individual rules for each
of the direct dependencies of the manuscript, which get assigned the
lowest possible priority in the main Snakefile. If one of these rules 
gets executed (when there are no other rules that can generate that
dependency), they simply throw an error informing the user about the problem.

"""
from showyourwork import exceptions


for dep in config["dependencies"][config["ms_tex"]]:
    rule:
        name:
            f"syw__fallback__{dep}"
        output:
            dep
        run:
            raise exceptions.MissingDependencyError(f"No rule to generate {dep}.")