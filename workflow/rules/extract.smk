rule extract:
    input:
        "src/data/results.tar.gz"
    output:
        "src/data/results_0.dat",
        "src/data/results_1.dat",
        "src/data/results_2.dat",
        "src/data/results_3.dat",
        "src/data/results_4.dat",
        "src/data/results_5.dat",
        "src/data/results_6.dat",
        "src/data/results_7.dat",
        "src/data/results_8.dat",
        "src/data/results_9.dat",
        
    shell:
        "tar -xzvf {input}"
