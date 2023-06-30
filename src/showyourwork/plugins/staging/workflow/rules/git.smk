from snakemake_staging.utils import package_data, working_directory

checkout = working_directory("staging", "git", "checkout", "{repo}", config=config)

rule staging__git__clone:
    output:
        directory(checkout)
    conda:
        package_data("workflow", "envs", "git.yml")
    params:
        url=lambda wildcards: get,
        branch="master"
    shell:
        "git clone {params.url} {output}"
