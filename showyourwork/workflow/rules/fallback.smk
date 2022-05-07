"""
Defines the ``syw__fallback__`` rules.

These are fallback rules for manuscript dependencies, which get invoked
if no other valid rule is found. Implementing these is necessary because 
of certain edge cases involving the behavior of Snakemake when dependencies 
aren't present.

Specifically, this Snakefile implements individual rules for each
of the direct dependencies of the manuscript, which get assigned the
lowest possible priority in the main Snakefile. If one of these rules 
gets executed (when there are no other rules that can generate that
dependency), the default behavior is to simply throw an error informing 
the user about the problem.

"""
from showyourwork import paths, exceptions


if config["fallback_rules"]:
    
    # Fail informatively
    for dep in config["dependencies"][config["ms_tex"]]:

        rule:
            name:
                f"syw__fallback__{dep}"
            output:
                dep
            run:
                raise exceptions.MissingDependencyError(
                    f"No rule to generate {dep}. This usually happens when "
                    "one or more of the dependencies of this file are missing "
                    "(and Snakemake doesn't know how to generate them). If "
                    "you're not sure why this error is being thrown, try "
                    "setting `fallback_rules: false` in the config file to get "
                    "a more informative error from Snakemake about possible "
                    "missing dependencies."
                )