<p align="center">
<a href="https://github.com/rodluger/showyourwork">
<img width = "450" src="https://raw.githubusercontent.com/rodluger/showyourwork/img/showyourwork.png" alt="showyourwork"/>
</a>
<br>
<br>
<a href="https://github.com/{{ GITHUB_SLUG }}/actions/workflows/showyourwork.yml">
<img src="https://github.com/{{ GITHUB_SLUG }}/actions/workflows/showyourwork.yml/badge.svg" alt="Article status"/>
</a>
<a href="https://{{ GITHUB_USER }}.github.io/{{ GITHUB_REPO }}/report.html">
<img src="https://img.shields.io/badge/workflow-report-blue.svg?style=flat" alt="View the workflow report"/>
</a>
<a href="https://github.com/{{ GITHUB_SLUG }}/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
</p>

## Your article is building

Your new repository is set up and ready to go. A GitHub Action is currently running that will generate your article PDF; you can check on the progress of the build [here](https://github.com/{{ GITHUB_SLUG }}/actions/workflows/showyourwork.yml). Once it's done, you can click on the badges at the top to take you to the workflow logs, the article graph, and the article PDF. The PDF is automatically updated every time you push changes to this repo; note that builds usually take a few minutes (or more, depending on what you're doing).

The first thing you might want to do is customize the `tex/ms.tex` file, which is currently just filled with placeholder text. You should also delete the current placeholder figure scripts in the `figures` directory, and add the scripts needed to build your own figures. If your workflow has external dependencies (which it most likely will), you must add them to the `environment.yml` file so `showyourwork` can build the paper from scratch. See [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#managing-environments) for details.

Finally, change or edit the `LICENSE` as needed and replace the text in this `README.md` with some useful information for the reader!

If you run into any trouble, please check out the [showyourwork documentation](https://showyourwork.readthedocs.io). If you think you've encountered a bug, please check out the [issues page](https://github.com/rodluger/showyourwork/issues) and raise a new one if needed.
