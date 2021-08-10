<p align="center">
<a href="https://github.com/rodluger/showyourwork">
<img width = "450" src="https://raw.githubusercontent.com/rodluger/showyourwork/main/showyourwork.png" alt="showyourwork"/>
</a>
<br>
<br>
{% if readme_message == "example" %}
An open source scientific article generated with <a href='https://github.com/rodluger/showyourwork'>showyourwork</a>.
{% elif readme_message == "test" %}
This is a test repository for <a href='https://github.com/rodluger/showyourwork'>showyourwork</a>. Not much to see here!
{% elif readme_message == "new" %}
Sit tight while your article builds. This may take up to five minutes, at which point you can refresh this page. You can also check the build progress by clicking on the Actions tab at the top. If there are any issues, please refer to the <a href='https://github.com/rodluger/showyourwork'>showyourwork documentation</a>.
{% endif %}
{% if not repo_active %}
<!--
{% endif %}
<p align="center">
<a href="https://github.com/{{ slug }}/actions/workflows/showyourwork.yml">
<img src="https://github.com/{{ slug }}/actions/workflows/showyourwork.yml/badge.svg" alt="Article status"/>
</a>
<a href="https://github.com/{{ slug }}/raw/main-pdf/dag.pdf">
<img src="https://img.shields.io/badge/workflow-graph-blue.svg?style=flat" alt="View the workflow DAG"/>
</a>
<a href="https://github.com/{{ slug }}/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
</p>
{% if not repo_active %}
-->
{% endif %}
</p>
