rule fibonacci_data:
    output:
        "data/fibonacci.dat"
    params:
        scripts=["figures/fibonacci.py"]
    shell:
        "curl https://zenodo.org/record/5187276/files/fibonacci.dat --output {output[0]}"