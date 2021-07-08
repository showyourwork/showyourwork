# article

If you don't have `mamba` installed:

```
conda install -c conda-forge mamba
```

Create and activate a new environment:

```
mamba create -p ./.env
mamba env update -p ./.env -f environment.yml
conda activate ./.env
```