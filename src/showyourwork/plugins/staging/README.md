# "Staging" for Snakemake

This package provides a mechanism for
[Snakemake](https://snakemake.readthedocs.io) workflows to explicitly "stage
out" the output files from certain rules to a public repository like
[Zenodo](https://zenodo.org) to allow faster re-execution of the workflow, using
these previously generated artifacts. This can be especially useful for
workflows with computationally expensive rules that don't need to be frequently
re-run.

`snakemake-staging` is a spin-off of the
[`showyourwork`](https://github.com/showyourwork/showyourwork) project, which
provides a "caching" framework for Snakemake workflows, to transparently avoid
re-execution of rules that have been cached to [Zenodo](https://zenodo.org). The
implementation of this logic in `showyourwork` is, however, somewhat fragile and
unpredictable. In `snakemake-staging`, we take a more explicit approach, where
"staged" rules are always either explicitly executed or restored.

## Installation

To use `snakemake-staging` in your workflow, you can install it using `pip`
(it's probably best to set up your Snakemake installation [following the
Snakemake
docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
first):

```bash
python -m pip install snakemake-staging
```

## Quickstart

### The `Snakefile`

While testing, it's probably best to use the [Zenodo
Sandbox](https://sandbox.zenodo.org), rather than the main site, since any
archive published to Zenodo is permanent. To use the sandbox, you'll need a
personal access token stored in the `SANDBOX_TOKEN` environment variable. You
can generate a new token
[here](https://sandbox.zenodo.org/account/settings/applications/).

Once you've added this token to your environment, you can edit the Snakefile for
your workflow to use `snakemake-staging` as follows. First, towards the top of
your Snakefile, add:

```python
from showyourwork.plugins import staging

stage = staging.ZenodoStage(
    "zenodo-stage",
    config.get("restore", False)
)
```

to create a new stage called `zenodo-stage`. Note that here we're extracting a
`restore` flag from the Snakemake config, which will be used to determine
whether to restore files for the stage. This means that you can control the
behavior of this stage from the command line. By passing `--config restore=True`
to the `snakemake` command line interface, all files staged out by the
`zenodo-stage` stage will be restored from the archive rather than generated.

Then, to stage out a rule, you can apply the stage as follows:

```python
rule expensive:
    input:
        ...
    output:
        stage(
            "path/to/output1.txt",
            "path/to/output2.txt",
        )
    shell:
        ...
```

Finally, _after defining all the rules that you want to stage out_, you must
add the following `include` which defines all the staging rules:

```python
include: staging.snakefile()
```

At this point, here's the full `Snakefile`:

<details>
<summary>Full Snakefile</summary>

```python
from showyourwork.plugins import staging

stage = staging.ZenodoStage(
    "zenodo-stage",
    config.get("restore", False)
)

rule expensive:
    input:
        ...
    output:
        stage(
            "path/to/output1.txt",
            "path/to/output2.txt",
        )
    shell:
        ...

include: staging.snakefile()
```

</details>

### Usage

With the `Snakefile` defined in the previous section, you can now run your
workflow in 3 ways:

1. **Normal execution**: If you run something like
   `snakemake path/to/output1.txt` (where I have omitted the usual `--cores`
   and `--conda` arguments) will execute the workflow as normal, without
   staging out any files.

2. **Stage upload**: If you instead have Snakemake target the `staging__upload`
   rule, the `expensive` rule will be executed, and the outputs will be uploaded
   to Zenodo, saving the record information to `zenodo-stage.zenodo.json` (this
   filename can be changed by passing the `info_file` argument to the
   `ZenodoStage` constructor).

3. **Stage restore**: Finally, after these outputs have been uploaded to Zenodo,
   you can call Snakemake `--config restore=True` to disable the `expensive`
   rule, and force the outputs to be restored from Zenodo.
