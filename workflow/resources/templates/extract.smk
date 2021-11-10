rule {{ rulename }}:
    message:
        "Extracting tarball `{{ input }}`..."
    input:
        "{{ input }}"
    output:
        {% for file in contents -%}
        "{{ file }}",
        {% endfor %}
    conda:
        posix(abspaths.user / "environment.yml")
    shell:
        "tar -xzvf {input}"