from pathlib import Path
import shutil


rule remove_zenodo:
    """
    On CI, we delete all datasets at the end so we don't cache them;
    this ensures we're always generating the figures based on the
    latest version of the deposit.
    
    """
    run:
    
        # Direct dependencies
        for file in files.zenodo_files_auto + files.zenodo_files_manual:
            if Path(file).exists():
                Path(file).unlink()

        # Contents of tarball dependencies
        for contents in zenodo.deposit_contents.values():
            for file in contents:
                if Path(file).exists():
                    if Path(file).is_dir():
                        shutil.rmtree(file)
                    else:
                        Path(file).unlink()