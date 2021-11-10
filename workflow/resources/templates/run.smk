rule {{ rulename }}:
    message:
        "Generating dataset(s) from file `{{ input }}`..."
    input:
        "{{ input }}"
    output:
        {% for file in contents -%}
        "{{ file }}",
        {% endfor %}
    conda:
        posix(abspaths.user / "environment.yml")
    shell:
        "{{ shell_cmd }}"