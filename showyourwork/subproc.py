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
    args, shell=False, cwd=None, secrets=[], callback=process_run_result
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

    """
    # Run the command and capture all output
    result = subprocess.run(
        args,
        shell=shell,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    # Parse the output
    stdout = result.stdout.decode()
    stderr = result.stderr.decode()
    code = result.returncode

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
    except:
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
        raise exceptions.RequestError(
            status=data["status"], message=data["message"]
        )
    else:
        return data
