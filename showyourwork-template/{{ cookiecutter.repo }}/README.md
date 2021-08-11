<p align="center">
<a href="https://github.com/rodluger/showyourwork">
<img width = "450" src="https://raw.githubusercontent.com/rodluger/showyourwork/main/showyourwork.png" alt="showyourwork"/>
</a>
<br>
<br>
{% if cookiecutter.repo_active == "y" %}
An open source scientific article generated with <a href='https://github.com/rodluger/showyourwork'>showyourwork</a>.
{% else %}
Sit tight while your article builds. This may take up to five minutes, at which point you can refresh this page. You can also check the build progress by clicking on the Actions tab at the top. If there are any issues, please refer to the <a href='https://github.com/rodluger/showyourwork'>showyourwork documentation</a>.
<!--
{% endif %}
<p align="center">
<a href="https://github.com/{{ cookiecutter.slug }}/actions/workflows/showyourwork.yml">
<img src="https://github.com/{{ cookiecutter.slug }}/actions/workflows/showyourwork.yml/badge.svg" alt="Article status"/>
</a>
<a href="https://github.com/{{ cookiecutter.slug }}/raw/main-pdf/dag.pdf">
<img src="https://img.shields.io/badge/workflow-graph-blue.svg?style=flat" alt="View the workflow DAG"/>
</a>
<a href="https://github.com/{{ cookiecutter.slug }}/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
</p>
{% if cookiecutter.repo_active == "n" %}

## You're all set!

Your new repository is set up and ready to go. The badges at the top will take you to the workflow logs, the article graph, and the compiled article PDF. The PDF is automatically updated every time you push changes to this repo; note that builds usually take a few minutes (or more, depending on what you're doing).

The first thing you might want to do is customize the `tex/ms.tex` file, which is currently just filled with placeholder text. You should also delete the current placeholder figure script in the `figures` directory, and add the scripts needed to build your own figures. If your workflow has external dependencies (which it most likely will), you must add them to the `environment.yml` file so `showyourwork` can build the paper from scratch. See [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#managing-environments) for details.

Finally, change or edit the `LICENSE` as needed and replace the text in this `README.md` with some useful informatin for the reader!

If you run into any trouble, please check out the [showyourwork documentation](https://github.com/rodluger/showyourwork). If you think you've encountered a bug, please check out the [issues page](https://github.com/rodluger/showyourwork/issues) and raise a new one if needed.
-->
{% endif %}
</p>