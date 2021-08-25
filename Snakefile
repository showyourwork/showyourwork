module showyourwork:
    snakefile:
        ".showyourwork/Snakefile"
    config:
        config


use rule * from showyourwork


rule fibonacci_data:
    output:
        "data/fibonacci.dat",
    shell:
        "curl https://zenodo.org/record/5187276/files/fibonacci.dat --output {output[0]}"


use rule figure from showyourwork as fibonacci_figure with:
    input:
        "figures/fibonacci.py",
        "data/fibonacci.dat",
        "environment.yml",
    output:
        "figures/fibonacci.pdf",


use rule figure from showyourwork as inline_figure with:
    input:
        "figures/inline.py",
        "environment.yml",
    output:
        "figures/inline.pdf",
