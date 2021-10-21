<p align="center">
<a href="https://github.com/rodluger/showyourwork">
<img width = "450" src="https://raw.githubusercontent.com/rodluger/showyourwork/img/showyourwork.png" alt="showyourwork"/>
<br>
<br>
<a href="https://github.com/rodluger/showyourwork/releases/tag/v0.1.13">
    <img src="https://img.shields.io/static/v1?label=version&message=0.1.13&color=blue"/>
</a>
<a href="https://showyourwork.readthedocs.io/en/v0.1.13">
    <img src="https://img.shields.io/static/v1?label=read&message=the%20docs&color=blue"/>
</a>
<a href="https://github.com/rodluger/showyourwork-template/generate">
    <img src="https://img.shields.io/static/v1?label=create&message=new%20repo&color=brightgreen"/>
</a>
</p>

<h2>Overview</h2>

<p align="justify">
This repository is intended to help authors publish the code that generated the figures and results in a scientific article. It ensures that the compiled article PDF is always in sync with all of the code used to generate it. It does this automatically—and seamlessly—with the help of the <a href="https://snakemake.readthedocs.io">Snakemake</a> workflow management system, the <a href="https://tectonic-typesetting.github.io">tectonic</a> typesetting engine, and <a href="https://github.com/features/actions">Github Actions CI</a>. The basic philosophy behind <code>showyourwork</code> is this: scientific papers should exist as GitHub repositories comprised of LaTeX files, figure scripts, rules to access datasets, a platform/environment specification, <i>and nothing else</i>. Anyone should be able to re-generate the article PDF from scratch at the click of a button.
</p>

<p align="justify">
Click <a href="https://github.com/rodluger/showyourwork-template/generate">here</a> to get started with a fresh article repository based on <code>showyourwork</code>. Once your repo is created, a GitHub Action will automatically run to finish setting it up. Refresh the page after a few minutes to view the new <code>README.md</code> with instructions on how to get started.
</p>

<h2>Test Suite</h2>

<p align="justify">
The <a href="https://github.com/rodluger/showyourwork-template">showyourwork-template</a> repository contains several branches, each containing a workflow template that demonstrates a particular use case or customization of <code>showyourwork</code>.
Upon every release of <code>showyourwork</code> (triggered via a GitHub Actions <a href="https://github.com/rodluger/showyourwork/actions/workflows/release.yml">workflow_dispatch</a> event), these templates are updated and instantiated
on the <a href="https://github.com/rodluger/showyourwork-example">showyourwork-example</a> repository. The following table
shows the build status of all of these examples.
</p>

<table>
  <tr>
    <th>branch</th>
    <th>build status</th>
    <th>output</th>
  </tr>
  
  <tr>
    <td>
        <a href="https://github.com/rodluger/showyourwork-example/tree/custom-figure-rule">custom-figure-rule</a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml?query=branch%3Acustom-figure-rule">
        <img src="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml/badge.svg?branch=custom-figure-rule" alt="test status"/>
    </a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/raw/custom-figure-rule-pdf/ms.pdf">
    <img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
    </a>
    </td>
  </tr>
  
  <tr>
    <td>
        <a href="https://github.com/rodluger/showyourwork-example/tree/expensive-figure">expensive-figure</a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml?query=branch%3Aexpensive-figure">
        <img src="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml/badge.svg?branch=expensive-figure" alt="test status"/>
    </a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/raw/expensive-figure-pdf/ms.pdf">
    <img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
    </a>
    </td>
  </tr>
  
  <tr>
    <td>
        <a href="https://github.com/rodluger/showyourwork-example/tree/figure-dataset">figure-dataset</a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml?query=branch%3Afigure-dataset">
        <img src="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml/badge.svg?branch=figure-dataset" alt="test status"/>
    </a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/raw/figure-dataset-pdf/ms.pdf">
    <img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
    </a>
    </td>
  </tr>
  
  <tr>
    <td>
        <a href="https://github.com/rodluger/showyourwork-example/tree/figure-deps">figure-deps</a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml?query=branch%3Afigure-deps">
        <img src="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml/badge.svg?branch=figure-deps" alt="test status"/>
    </a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/raw/figure-deps-pdf/ms.pdf">
    <img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
    </a>
    </td>
  </tr>
  
  <tr>
    <td>
        <a href="https://github.com/rodluger/showyourwork-example/tree/main">main</a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml?query=branch%3Amain">
        <img src="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml/badge.svg?branch=main" alt="test status"/>
    </a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/raw/main-pdf/ms.pdf">
    <img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
    </a>
    </td>
  </tr>
  
  <tr>
    <td>
        <a href="https://github.com/rodluger/showyourwork-example/tree/multi-panel-figure">multi-panel-figure</a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml?query=branch%3Amulti-panel-figure">
        <img src="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml/badge.svg?branch=multi-panel-figure" alt="test status"/>
    </a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/raw/multi-panel-figure-pdf/ms.pdf">
    <img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
    </a>
    </td>
  </tr>
  
  <tr>
    <td>
        <a href="https://github.com/rodluger/showyourwork-example/tree/multiple-figures">multiple-figures</a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml?query=branch%3Amultiple-figures">
        <img src="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml/badge.svg?branch=multiple-figures" alt="test status"/>
    </a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/raw/multiple-figures-pdf/ms.pdf">
    <img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
    </a>
    </td>
  </tr>
  
  <tr>
    <td>
        <a href="https://github.com/rodluger/showyourwork-example/tree/simple-figure">simple-figure</a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml?query=branch%3Asimple-figure">
        <img src="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml/badge.svg?branch=simple-figure" alt="test status"/>
    </a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/raw/simple-figure-pdf/ms.pdf">
    <img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
    </a>
    </td>
  </tr>
  
  <tr>
    <td>
        <a href="https://github.com/rodluger/showyourwork-example/tree/static-figure">static-figure</a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml?query=branch%3Astatic-figure">
        <img src="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml/badge.svg?branch=static-figure" alt="test status"/>
    </a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/raw/static-figure-pdf/ms.pdf">
    <img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
    </a>
    </td>
  </tr>
  
</table>
