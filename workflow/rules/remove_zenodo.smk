from pathlib import Path


rule remove_zenodo:
    """
    On CI, we delete all datasets at the end so we don't cache them;
    this ensures we're always generating the figures based on the
    latest version of the deposit.
    
    """
    run:
        for file in files.zenodo_files_auto:
            if Path(file).exists():
                Path(file).unlink()