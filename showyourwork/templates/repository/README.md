<p align="center">
  <a href="https://github.com/rodluger/showyourwork">
      <img width = "450" src="https://raw.githubusercontent.com/rodluger/showyourwork/main/showyourwork.png" alt="showyourwork"/>
  </a>
  <br>
  <br>
  {{ readme_message }}
  <p align="center">
    {% if repo_active %}
    <a href="https://github.com/{{ slug }}/actions/workflows/showyourwork.yml">
      <img src="https://github.com/{{ slug }}/actions/workflows/showyourwork.yml/badge.svg" alt="Article status"/>
    </a>
    <a href="https://github.com/{{ slug }}/raw/main-pdf/dag.pdf">
      <img src="https://img.shields.io/badge/workflow-graph-blue.svg?style=flat" alt="View the workflow DAG"/>
    </a>
    <a href="https://github.com/{{ slug }}/raw/main-pdf/ms.pdf">
      <img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
    </a>
    {% else %}
    <a href="#">
      <img src="https://img.shields.io/badge/build-in%20progress-yellow.svg?style=flat" alt="Article status"/>
    </a>
    {% endif %}
</p>
