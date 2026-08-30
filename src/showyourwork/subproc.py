import os
import signal
import subprocess


def process_run_result(code, stdout, stderr):
    """
    Default callback function for ``get_stdout``.

    """
    from . import exceptions, logging

    # Log the output
    logger = logging.get_logger()
    if stdout:
        logger.debug(stdout)

    # Raise the exception
    if code != 0:
        raise exceptions.CalledProcessError(stderr)

    return stdout


def get_stdout(
    args,
    shell=False,
    cwd=None,
    secrets=(),
    callback=process_run_result,
    env=None,
    capture_output=True,
    timeout=None,
):
    """
    A thin wrapper around ``subprocess.run`` that hides secrets and decodes
    ``stdout`` and ``stderr`` output into ``utf-8``.

    Args:
        args (list or str): Arguments passed to ``subprocess.run``
        shell (bool, optional): Passed directly to ``subprocess.run``
        cwd (str, optional): Directory to run the command in, if different
            from current working directory.
        secrets (list, optional): Secrets to be masked in the output.
        callback (callable, optional): Callback to process the result.
        env (dict): Extra environment variables to be passed to the ``subprocess.run``
        capture_output (bool, optional): If ``True`` (default), capture ``stdout`` and
            ``stderr`` and pass decoded strings to ``callback``. If ``False``, stream
            command output directly to the terminal.
        timeout (float, optional): Maximum number of seconds to let the command
            run before it is killed. ``None`` (default) disables the timeout,
            preserving the previous blocking behavior. Without this, a stalled
            child process (e.g. a build hanging on a network call several
            subprocess layers down) can block indefinitely until an external
            watchdog such as GitHub Actions' 6-hour job cap kills the job.

    """
    #  Update the environment variables if passed
    subprocess_env = os.environ.copy()
    if env is not None:
        subprocess_env.update(env)

    # Run the command in its own process group so that, on timeout, we can
    # kill the whole subtree it spawned (e.g. `sh -c "showyourwork build"`
    # -> snakemake -> a rule's python subprocess) rather than just the
    # immediate child. subprocess.run's built-in timeout handling only kills
    # the direct child object it manages, which with shell=True leaves
    # grandchildren running as orphans -- exactly what we saw GitHub Actions'
    # post-job cleanup having to mop up after the 6-hour job cap killed the
    # job. Using Popen directly lets us os.killpg() the whole group ourselves.
    popen_kwargs = {}
    if os.name != "nt":
        popen_kwargs["start_new_session"] = True
    else:
        popen_kwargs["creationflags"] = subprocess.CREATE_NEW_PROCESS_GROUP

    process = subprocess.Popen(
        args,
        shell=shell,
        cwd=cwd,
        env=subprocess_env,
        stdout=subprocess.PIPE if capture_output else None,
        stderr=subprocess.PIPE if capture_output else None,
        **popen_kwargs,
    )
    try:
        stdout, stderr = process.communicate(timeout=timeout)
        code = process.returncode
    except subprocess.TimeoutExpired:
        if os.name != "nt":
            try:
                os.killpg(os.getpgid(process.pid), signal.SIGKILL)
            except ProcessLookupError:
                pass
        else:
            process.kill()
        # Drain any output the process had already produced before it was
        # killed, so it isn't lost from the error message below.
        stdout, stderr = process.communicate()
        code = 124  # conventional exit code for a timed-out command

    # Parse the output
    stdout = stdout.decode() if isinstance(stdout, bytes) else stdout
    stderr = stderr.decode() if isinstance(stderr, bytes) else stderr
    stdout = stdout or ""
    stderr = stderr or ""
    if code == 124:
        stderr += f"\n[Command timed out after {timeout} seconds: {args!r}]"

    # Hide secrets from the command output
    for secret in secrets:
        stdout = stdout.replace(secret, "*****")
        stderr = stderr.replace(secret, "*****")

    # Callback
    return callback(code, stdout, stderr)


def parse_request(r):
    """
    Parse a requests return object ``r`` and raise a custom exception
    for a >200-level status code.

    """
    from . import exceptions

    # Try to get the data
    try:
        data = r.json()
    except Exception:
        if len(r.text) == 0:
            # We're good; there's just no data
            data = {}
        else:
            # Something went wrong
            data = {"message": r.text}

    # Parse the status code
    if r.status_code > 204:
        data["message"] = data.get("message", "")
        data["status"] = data.get("status", "")
        for error in data.get("errors", []):
            data["message"] += " " + error.get("message", "")
        raise exceptions.RequestError(status=data["status"], message=data["message"])
    else:
        return data
