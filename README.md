<p align="center">
  <img width = "450" src=".showyourwork/tex/showyourwork.png"/>
  <br>
  <br>
  <a href="https://github.com/rodluger/showyourwork/generate">
    <img src="https://img.shields.io/badge/Use%20this%20template-2EA44F.svg?style=flat" alt="Use this template"/>
  </a>
  <a href="https://github.com/rodluger/showyourwork/blob/main-pdf/ms.pdf">
    <img src="https://img.shields.io/badge/Read%20sample%20article-blue.svg?style=flat" alt="Use this template"/>
  </a>
</p>


### Overview

`showyourwork` is a [repository template](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/creating-a-repository-on-github/creating-a-repository-from-a-template) aimed at helping you write fully reproducible scientific articles in LaTeX. 

### Setup

Click [here](https://github.com/rodluger/showyourwork/generate) to create your own repository based on `showyourwork` (you can name it whatever you'd like). Once you clone it locally and `cd` into your repo, run

```
conda install -c conda-forge mamba
mamba create -p ./envs
mamba env update -p ./envs -f environment.yml
conda activate ./envs
```

to set up the default conda environment for the project (if you don't have conda, install it [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)).

### Build locally

Compile the article with `snakemake` (installed automatically in the step above):

```
snakemake -c1
```

This will generate `ms.pdf` at the root of the repository. The article is built from a manuscript file [tex/ms.tex](https://github.com/rodluger/showyourwork/blob/main/tex/ms.tex), an optional bibliography file [tex/bib.bib](https://github.com/rodluger/showyourwork/blob/main/tex/bib.bib), and any number of Python figure scripts in the [figures](https://github.com/rodluger/showyourwork/tree/main/figures) directory and any number of Python equation test scripts in the [tests](https://github.com/rodluger/showyourwork/tree/main/tests) directory.

### Build in the cloud

This repository is set up to automatically compile the paper upon every commit with [GitHub Actions](https://github.com/rodluger/showyourwork/actions). The compiled [ms.pdf](https://github.com/rodluger/showyourwork/blob/main-pdf/ms.pdf) is automatically made available on the [main-pdf](https://github.com/rodluger/showyourwork/tree/main-pdf) branch.