from ... import paths

# A wrapper around snakemake; has the same effect as calling `snakemake` from
# the command line, but allows us to apply monkeypatches before the workflow starts.
SNAKEMAKE = f"python {paths.showyourwork().module / 'cli' / 'patches' / 'launch_snakemake.py'}"
