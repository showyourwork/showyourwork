import subprocess
import os
from pathlib import Path


def process_run_result(code, stdout, stderr):
    from . import logging, exceptions

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
    Parse a requests return object and raise a custom exception
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
            data["message"] += " " + error["message"]
        raise exceptions.RequestError(
            status=data["status"], message=data["message"]
        )
    else:
        return data