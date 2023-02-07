from showyourwork import paths

rule syw__default_copy_manuscript:
    input:
        paths.repo(config).manuscript
    output:
        paths.work(config).manuscript
    run:
        import shutil
        shutil.copyfile(input[0], output[0])
