"""
This script applies various monkeypatches to Snakemake and then invokes
it as if it were run from the command line.

"""


def patch_check_and_touch_output():
    """
    Fixes a backwards-compatibility issue first identified here:

        https://github.com/showyourwork/showyourwork/pull/200#issuecomment-1259757992

    Showyourwork monkeypatches the `wait_for_files` method of the Snakemake
    DAG in order to display a custom error message. Recently, however, we
    upgraded the pinned Snakemake version from 6.15.5 to 7.14.0. In this newer
    version of Snakemake, one of the keyword arguments to `wait_for_files`
    changed names from `ignore_pipe` to `ignore_pipe_or_service`. We've upgraded
    the monkeypatch in `showyourwork/patches.py`, but there's still a nasty
    backwards compatibility issue: if a workflow specifying a version of
    showyourwork <= 0.3.1 is executed with showyourwork > 0.3.1, the workflow
    will apply the _wrong_ monkeypatch and a `TypeError` will be thrown saying
    either `ignore_pipe` or `ignore_pipe_or_service` is not a valid keyword
    to `wait_for_files`.

    The solution is a _second_ monkeypatch on the function `check_and_touch_output`,
    which is currently (as of 7.14.0) the only instance in the codebase where
    `wait_for_files` is explicitly called with the `ignore_pipe_or_service`
    keyword. Here we simply patch it with a `try-except` block to try both
    keyword names if needed.

    """

    def patch(method):
        def patched_method(*args, **kwargs):
            def wait_for_files(*args, **kwargs):
                ignore_pipe = kwargs.pop("ignore_pipe", False) or kwargs.pop(
                    "ignore_pipe_or_service", False
                )
                try:
                    return snakemake.io.wait_for_files(
                        *args, ignore_pipe=ignore_pipe, **kwargs
                    )
                except TypeError:
                    return snakemake.io.wait_for_files(
                        *args, ignore_pipe_or_service=ignore_pipe, **kwargs
                    )

            snakemake.dag.wait_for_files = wait_for_files
            return method(*args, **kwargs)

        return patched_method

    snakemake.dag.DAG.check_and_touch_output = patch(
        snakemake.dag.DAG.check_and_touch_output
    )


if __name__ == "__main__":

    import snakemake

    # Patch the `check_and_touch_output` method of the DAG
    patch_check_and_touch_output()

    # Launch snakemake. Any command line arguments are automatically forwarded.
    snakemake.main()
