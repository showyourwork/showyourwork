<p align="center">
  <img width = "450" src=".showyourwork/tex/showyourwork.png"/>
</p>


### Overview

`showyourwork` is a repository template aimed at helping you write fully reproducible scientific articles in LaTeX.

### Requirements

`showyourwork` requires both `conda` and `mamba`.

If you don't have `conda` installed, see [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for instructions. And if you don't have `mamba`, you can easily install it with `conda`:

```
conda install -c conda-forge mamba
```

### Setup

This is a [repository template](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/creating-a-repository-on-github/creating-a-repository-from-a-template). Click [here](https://github.com/rodluger/showyourwork/generate) to create your own repository based on `showyourwork`, name it whatever you'd like, and clone it locally. Then, activate the `mamba` environment within the repository as follows:

```
mamba create -p ./.env
mamba env update -p ./.env -f environment.yml
conda activate ./.env
```

### Build

Compile the paper with `snakemake` (installed automatically in the step above):

```
snakemake -c1
```

This will generate `ms.pdf` at the root of the repository.

### Continuous integration

This repository is set up to automatically compile the paper upon every commit with [GitHub Actions](https://github.com/rodluger/showyourwork/actions). The compiled [ms.pdf](https://github.com/rodluger/showyourwork/blob/main-pdf/ms.pdf) is available on the `main-pdf` branch.