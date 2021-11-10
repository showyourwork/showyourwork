rule {{ rulename }}:
    message:
        "Generating tarball `{{ output }}`..."
    input:
        {% for file in contents -%}
        "{{ file }}",
        {% endfor %}
    output:
        "{{ output }}"
    conda:
        posix(abspaths.user / "environment.yml")
    shell:
        "tar -czvf {output} {input}"