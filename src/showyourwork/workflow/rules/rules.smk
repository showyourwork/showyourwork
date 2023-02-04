rule _syw_default_copy_manuscript:
    input:
        paths.repo(config).manuscript
    output:
        paths.work(config).manuscript
    shell:
        """
        cp "{input}" "{output}"
        """
