from snakemake_staging import stages, utils

for name, stage in stages.STAGES.items():
    if not isinstance(stage, stages.NoOpStage):
        continue

    # Rules for restoring or snapshotting the staging directory based on the
    # restore configuration
    if stage.restore:
        rule:
            name:
                utils.rule_name("noop", name)
            message:
                f"Restoring staging directory for '{name}'"
            output:
                [stage.directory / f for f in stage.files.keys()]
            run:
                # We check to make sure that all the files were restored
                for file in output:
                    if not file.exists():
                        raise RuntimeError(
                            f"File '{file}' was not sucessfully restored by "
                            f"stage '{name}'"
                        )

    else:
        rule:
            name:
                utils.rule_name("noop", name)
            message:
                f"Snapshotting stage '{name}'"
            input:
                [stage.directory / f for f in stage.files.keys()]
            output:
                touch(stage.upload_flag_file)

    # Rules for copying files to and from the staging directory based on the
    # restore configuration
    for staged_filename, filename in stage.files.items():
        if stage.restore:
            rule:
                name:
                    utils.rule_name("noop", name, "copy", path=filename)
                message:
                    f"Copying file '{filename}' from stage '{name}'"
                input:
                    stage.directory / staged_filename
                output:
                    filename
                run:
                    utils.copy_file_or_directory(input[0], output[0])

        else:
            rule:
                name:
                    utils.rule_name("noop", name, "copy", path=filename)
                message:
                    f"Copying file '{filename}' to stage '{name}'"
                input:
                    filename
                output:
                    stage.directory / staged_filename
                run:
                    utils.copy_file_or_directory(input[0], output[0])
