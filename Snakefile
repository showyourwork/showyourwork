# Import the showyourwork module
module showyourwork:
    snakefile:
        ".showyourwork/Snakefile"
    config:
        config


# Use all default rules
use rule * from showyourwork


# Subclass the `figure` rule to specify that `figures/koch.py`
# generates *two* figures
use rule figure from showyourwork as koch_figure with:
    input:
        "figures/koch.py",
        "environment.yml"
    output:
        report("figures/koch1.pdf", category="Figure"),
        report("figures/koch2.pdf", category="Figure")


# Custom rule to download a dataset
rule fibonacci_data:
    output:
        report("data/fibonacci.dat", category="Dataset")
    shell:
        "curl https://zenodo.org/record/5187276/files/fibonacci.dat --output {output[0]}"


# Subclass the `figure` rule to specify that `data/fibonacci.dat`
# is a dependency of `figures/fibonacci.pdf`
use rule figure from showyourwork as fibonacci_figure with:
    input:
        "figures/fibonacci.py",
        "data/fibonacci.dat",
        "environment.yml"
    output:
        report("figures/fibonacci.pdf", category="Figure")


# Subclass the `figure` rule to specify that the inline figure
# `figures/inline.pdf` is generated from the script `figures/inline.py`
use rule figure from showyourwork as inline_figure with:
    input:
        "figures/inline.py",
        "environment.yml"
    output:
        report("figures/inline.pdf", category="Figure")
