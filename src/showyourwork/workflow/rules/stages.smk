from showyourwork import stages

_restore = stages.get_stages_to_restore(config)
for stage, files in stages.STAGES.items():
    if name in _restore:
        rule:
            name:
                utils.rule_name("stages", "restore", stage)
            message:
                f"Restoring {len(files)} files from snapshot for stage '{stage}'"
            input:
                zenodo=stages.optionally_require_zenodo(config)
            output:
                files
            run:
                stages.restore_stage(config, stage)

    else:
        rule:
            name:
                utils.rule_name("stages", "snapshot", stage)
            message:
                f"Snapshotting {len(files)} files for stage '{stage}'"
            input:
                files,
                zenodo=stages.optionally_require_zenodo(config),
            output:
                touch(SYW__WORK_PATHS.flags(f"stage_{stage}.snapshot"))
            run:
                stages.snapshot_stage(config, stage)
