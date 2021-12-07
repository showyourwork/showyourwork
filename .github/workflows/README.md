<p align="center">
<a href="https://github.com/rodluger/showyourwork">
<img width = "450" src="https://raw.githubusercontent.com/rodluger/showyourwork/img/showyourwork.png" alt="showyourwork"/>
<br>
<br>
<a href="https://github.com/rodluger/showyourwork/releases/tag/v{{ version }}">
    <img src="https://img.shields.io/static/v1?label=version&message={{ version }}&color=blue"/>
</a>
<a href="https://showyourwork.readthedocs.io/en/v{{ version }}">
    <img src="https://img.shields.io/static/v1?label=read&message=the%20docs&color=blue"/>
</a>
<a href="https://github.com/rodluger/showyourwork-template/generate">
    <img src="https://img.shields.io/static/v1?label=create&message=new%20repo&color=brightgreen"/>
</a>
</p>

<h2>Overview</h2>
<p align="justify">
This repository is intended to help authors publish the code that generated the figures and results in a scientific article. It ensures that the compiled article PDF is always in sync with all of the code used to generate it, and that any user that clones the repository can reproduce the PDF <i>from scratch</i> by running
</p>

```
make
```

<p align="justify">
It does this automatically—and seamlessly—with the help of the <a href="https://snakemake.readthedocs.io">Snakemake</a> workflow management system, the <a href="https://www.anaconda.com/products/individual">conda</a> package manager, the <a href="https://tectonic-typesetting.github.io">tectonic</a> typesetting engine, the <a href="https://zenodo.org">Zenodo</a> data hosting service, and <a href="https://github.com/features/actions">Github Actions CI</a>. 
</p>

<h2>The showyourwork philosophy</h2>
<p align="justify">
Scientific papers should exist as GitHub repositories comprised of LaTeX files, figure scripts, rules to generate and access datasets, a platform/environment specification, <i>and nothing else</i>. Anyone should be able to re-generate the article PDF from scratch at the click of a button.
</p>

<h2>Getting started</h2>
<p align="justify">
Click <a href="https://github.com/rodluger/showyourwork-template/generate">here</a> to get started with a fresh article repository based on <code>showyourwork</code>. Name it whatever you'd like, and make sure to keep the <strong>Include all branches</strong> field unchecked.
Once your repo is created, a GitHub Action will automatically run to finish setting it up. Refresh the page after a few minutes to view the new <code>README.md</code> with instructions on how to get started. Please
check out <a href="https://showyourwork.readthedocs.io">the documentation</a> for more information, examples, tutorials, and FAQs.
</p>

<h2>Examples and test suite</h2>
<p align="justify">
The <a href="https://github.com/rodluger/showyourwork-template">showyourwork-template</a> repository contains several branches, each containing a workflow template that demonstrates a particular use case or customization of <code>showyourwork</code>.
Upon every release of <code>showyourwork</code> (triggered via a GitHub Actions <a href="https://github.com/rodluger/showyourwork/actions/workflows/release.yml">workflow_dispatch</a> event), these templates are updated and instantiated
on the <a href="https://github.com/rodluger/showyourwork-example">showyourwork-example</a> repository. The following table
shows the build status of all of these examples.
Please check out <a href="https://showyourwork.readthedocs.io/en/stable/custom/">the documentation</a> for more information about each one.
</p>

<table>
  <tr>
    <th>branch</th>
    <th>build status</th>
    <th>output</th>
  </tr>
  {% for branch in branches %}
  <tr>
    <td>
        <a href="https://github.com/rodluger/showyourwork-example/tree/{{ branch }}">{{ branch }}</a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml?query=branch%3A{{ branch }}">
        <img src="https://github.com/rodluger/showyourwork-example/actions/workflows/showyourwork.yml/badge.svg?branch={{ branch }}" alt="test status"/>
    </a>
    </td>
    <td>
    <a href="https://github.com/rodluger/showyourwork-example/raw/{{ branch }}-pdf/ms.pdf">
    <img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
    </a>
    </td>
  </tr>
  {% endfor %}
</table>

<h2>Attribution</h2>
<p align="justify">
We're working on writing up a citeable paper on <code>showyourwork</code>, but in the meantime, please
consider citing <a href="https://ui.adsabs.harvard.edu/abs/2021arXiv211006271L">this paper</a>,
in which we first introduced <code>showyourwork</code>:
</p>

```
@ARTICLE{Luger2021,
       author = {{Luger}, Rodrigo and {Bedell}, Megan and {Foreman-Mackey}, Daniel and {Crossfield}, Ian J.~M. and {Zhao}, Lily L. and {Hogg}, David W.},
        title = "{Mapping stellar surfaces III: An Efficient, Scalable, and Open-Source Doppler Imaging Model}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2021,
        month = oct,
          eid = {arXiv:2110.06271},
        pages = {arXiv:2110.06271},
archivePrefix = {arXiv},
       eprint = {2110.06271},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021arXiv211006271L},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

<p align="justify">
We'd also appreciate it if you'd let us know about your project so we can showcase it
as an example on the documentation; consider adding your project to our
<a href="https://github.com/rodluger/showyourwork/edit/main/docs/projects.json">list of curated repos</a>.
And if you're interested in contributing to <code>showyourwork</code>, please consider
opening a pull request to tackle one of the many <a href="https://github.com/rodluger/showyourwork/issues?q=is%3Aissue+is%3Aopen+">open issues</a> or open a new issue with requests/suggestions for new functionality.
</p>