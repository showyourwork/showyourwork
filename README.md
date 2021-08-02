<p align="center">
  <img width = "450" src="./showyourwork.png"/>
  <br>
  <br>
  <a href="https://github.com/rodluger/showyourwork/actions/workflows/test_dispatch.yml">
    <img src="https://github.com/rodluger/showyourwork/actions/workflows/test_dispatch.yml/badge.svg"/>
  </a>
  <a href="https://github.com/rodluger/showyourwork/actions/workflows/test_receive.yml">
    <img src="https://github.com/rodluger/showyourwork/actions/workflows/test_receive.yml/badge.svg"/>
  </a>
</p>

## Overview

This package is intended to help authors publish the code that generated the results (and in particular the figures) in a scientific article. It ensures that the compiled article PDF is always in sync with all of the code used to generate it. It does this automatically—and seamlessly—with the help of the [Snakemake](https://snakemake.readthedocs.io) workflow management system, the [AASTeX](https://journals.aas.org/aastexguide) markup package, the [tectonic](https://tectonic-typesetting.github.io) typesetting engine, and [Github Actions](https://github.com/features/actions) CI.

The basic philosophy behind `showyourwork` is this: scientific papers should exist as GitHub repositories comprised of LaTeX files, figure scripts, optional experiment datasets, _and nothing else_. Anyone should be able to re-generate the article PDF from scratch at the click of a button.

## Installation

`showyourwork` requires the [conda](https://www.anaconda.com/products/individual) package manager.

_NOTE: Currently only the development version is available._

```
conda install -c conda-forge -c bioconda mamba pip snakemake tectonic
python -m pip install git+https://github.com/rodluger/showyourwork.git@main
```

## Setup

To create a fresh article repository, run

```bash
showyourwork --new [version]
```

where `[version]` is optional and defaults to the latest version of this package. You will be prompted for some information about the repository, which will then be initialized in the current working directory.

## Repository structure

Your new repository should have the following structure:

```
my_repo
├── .github
│   └── workflows
│       └── showyourwork.yml
├── data
├── figures
│   ├── matplotlibrc
│   └── *.py
├── tex
│   ├── ms.tex
│   └── bib.bib
└── environment.yml
```

Let's discuss each file/directory one at a time:

- `.github/workflows/showyourwork.yml`: This is the GitHub Action workflow, which tells GitHub to recompile your article every time you commit and push changes to the remote repository. You shouldn't typically have to edit this file.

- `data`: This directory can contain optional data files required by your build. In general, these should be the results of experiments or observations—not something that can be generated programatically, like the result of a simulation (more on this below). Keep in mind that you shouldn't use `git` to track files larger than around `1 MB`; if your data set is very large, consider hosting it externally and using something like `wget` to access it, e.g., within a figure script.

- `figures`: This directory should contain all the scripts needed to generate the figures in your article. Currently, these must be Python (`*.py`) scripts, but support for other formats/languages is coming soon. We recommend one script per figure (or `\includegraphics` call) in the article, but this is not required (see below for details).

- `tex`: This is the directory containing your LaTeX manuscript (`ms.tex`) and optionally a bibliography file (`bib.bib`). You don't usually need anything else in here, unless your manuscript depends on custom stylesheets or classes not hosted on the traditional LaTeX channels. More details on the contents of `ms.tex` below.

- `environment.yml`: This contains the `conda` environment specs for building your article. All dependencies needed to build your article from scratch should be included here. For more details, see [the conda docs](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually).

## Building locally

To build your article locally, simply run

```bash
showyourwork
```

in the top level of your repository. This will generate the compiled PDF, `ms.pdf`, in the same directory. If this is the first time you're running it, `showyourwork` will build all your figure scripts before compiling the PDF. Once they are built, `showyourwork` caches the output to avoid re-running scripts that have not changed.

## Building remotely

The main point of `showyourwork` is to make your article open source, so the final step is to push it to GitHub. If you haven't already done so, run

```bash
git init
git add .
git commit -m "Initial commit"
```

from the top level to start tracking it with `git`. Create a [new repository](https://github.com/new) on GitHub, then push your changes:

```bash
git branch -M main
git remote add origin https://github.com/[user]/[repo].git
git push -u origin main
```

Upon every commit pushed to GitHub, `showyourwork` will **automatically build your article on GitHub Actions**. Upon completion, the PDF will be available both as a [workflow artifact](https://docs.github.com/en/actions/managing-workflow-runs/downloading-workflow-artifacts) and on a separate branch in your repo called `main-pdf` (or `[branch]-pdf` if you're on a different branch).

## Under construction...

This README is still under construction. Please stay tuned!
