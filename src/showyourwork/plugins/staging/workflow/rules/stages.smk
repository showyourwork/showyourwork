from pathlib import Path
from snakemake_staging import stages, utils
from snakemake_staging.config import _CONFIG

previously_included = set()
for name, stage in stages.STAGES.items():
    # Custom rules for uploading and downloading this type of stage
    snakefile = stage.snakefile()
    if snakefile not in previously_included:
        include: snakefile
        previously_included.add(snakefile)


# A rule to upload all stages
working_directory = Path(_CONFIG.get("working_directory", "staging"))
rule staging__upload:
    message:
        "Uploading all stages"
    input:
        [
            stage.upload_flag_file
            for stage in stages.STAGES.values()
            if not stage.restore
        ]
    output:
        touch(working_directory / "staging_upload.done")
