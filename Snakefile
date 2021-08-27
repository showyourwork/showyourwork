# Import the showyourwork module
module showyourwork:
    snakefile:
        ".showyourwork/Snakefile"
    config:
        config


# Use all default rules
use rule * from showyourwork


# Custom rule to download a dataset
rule fibonacci_data:
    output:
        report("src/figures/fibonacci.dat", category="Dataset")
    shell:
        "curl https://zenodo.org/record/5187276/files/fibonacci.dat --output {output[0]}"


# Subclass the `figure` rule to specify that `data/fibonacci.dat`
# is a dependency of `figures/fibonacci.pdf`
use rule figure from showyourwork as fibonacci_figure with:
    input:
        "src/figures/fibonacci.py",
        "src/figures/fibonacci.dat",
        "environment.yml"
    output:
        report("src/figures/fibonacci.pdf", category="Figure")


# Subclass the `figure` rule to specify that the inline figure
# `figures/inline.pdf` is generated from the script `figures/inline.py`
use rule figure from showyourwork as inline_figure with:
    input:
        "src/figures/inline.py",
        "environment.yml"
    output:
        report("src/figures/inline.pdf", category="Figure")
