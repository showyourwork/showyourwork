<!-- SHOWYOURWORK_TEMPLATE -->
<p align="center">
<a href="https://github.com/rodluger/showyourwork">
<img width = "450" src="https://raw.githubusercontent.com/rodluger/showyourwork/img/showyourwork.png" alt="showyourwork"/>
<br>
<br>
<a href="#">
    <img src="https://img.shields.io/static/v1?label=read&message=the%20docs&color=blue"/>
</a>
<a href="#">
    <img src="https://img.shields.io/static/v1?label=create&message=new%20repo&color=brightgreen"/>
</a>
</a>
</p>
<br>

`showyourwork` is intended to help authors publish the code that generated the results (and in particular the figures) in a scientific article. It ensures that the compiled article PDF is always in sync with all of the code used to generate it. It does this automatically—and seamlessly—with the help of the [Snakemake](https://snakemake.readthedocs.io) workflow management system, the [tectonic](https://tectonic-typesetting.github.io) typesetting engine, and [Github Actions](https://github.com/features/actions) CI. The basic philosophy behind `showyourwork` is this: scientific papers should exist as GitHub repositories comprised of LaTeX files, figure scripts, rules to access datasets, a platform/environment specification, _and nothing else_. Anyone should be able to re-generate the article PDF from scratch at the click of a button.
