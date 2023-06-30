from snakemake_staging import stages, utils, zenodo

for name, stage in stages.STAGES.items():
    if not isinstance(stage, zenodo.ZenodoStage):
        continue

    # Rules for restoring or snapshotting the staging directory based on the
    # restore configuration
    if stage.restore:
        for file in stage.files.values():
            rule:
                name:
                    utils.rule_name("zenodo", name, "download", path=file)
                message:
                    f"Restoring file '{file}' for stage '{name}'"
                input:
                    stage.info_file
                output:
                    file
                run:
                    stage.download_file(input[0], file)

    else:
        rule:
            name:
                utils.rule_name("zenodo", name, "draft")
            message:
                f"Creating draft for stage '{name}'"
            output:
                stage.draft_info_file
            run:
                stage.create_draft(output[0])

        for file in stage.files.values():
            rule:
                name:
                    utils.rule_name("zenodo", name, "upload", path=file)
                message:
                    f"Uploading file '{file}' for stage '{name}'"
                input:
                    stage.draft_info_file,
                    file
                output:
                    stage.upload_info_file(file)
                run:
                    stage.upload_file(input[0], input[1], output[0])

        rule:
            name:
                utils.rule_name("zenodo", name, "publish")
            message:
                f"Publishing stage '{name}'"
            input:
                stage.draft_info_file,
                [stage.upload_info_file(file) for file in stage.files.values()]
            output:
                stage.info_file,
                touch(stage.upload_flag_file)
            run:
                stage.publish_draft(input[0], output[0])
