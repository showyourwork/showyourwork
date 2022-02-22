from . import logging, exceptions
import subprocess


def process_run_result(code, stdout, stderr):
    # Log the output
    logger = logging.get_logger()
    if stdout:
        logger.debug(stdout)

    # Raise the exception with no traceback to hide
    # the invocation (which may contain secrets)
    if code != 0:
        with exceptions.no_traceback():
            raise exceptions.CalledProcessError(stderr)


def run(args, cwd=None, secrets=[], callback=process_run_result):
    # Run the command and capture all output
    result = subprocess.run(
        args,
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