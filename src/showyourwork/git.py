import subprocess
from functools import lru_cache
from typing import Any, Iterable, Optional, Union


def git(
    args: Union[Iterable[str], str], **kwargs: Any
) -> subprocess.CompletedProcess[str]:
    kwargs["check"] = kwargs.get("check", True)
    kwargs["text"] = kwargs.get("text", True)
    kwargs["capture_input"] = kwargs.get("capture_input", True)
    if isinstance(args, str):
        args = f"git {args}"
        kwargs["shell"] = kwargs.get("shell", True)
    else:
        args = ["git", *args]
        kwargs["shell"] = kwargs.get("shell", False)
    return subprocess.run(args, **kwargs)


@lru_cache
def get_repo_branch() -> str:
    return git(["rev-parse", "--abbrev-ref", "HEAD"]).stdout.strip()


@lru_cache
def get_repo_sha() -> str:
    return git(["rev-parse", "HEAD"]).stdout.strip()


@lru_cache
def get_repo_tag() -> Optional[str]:
    r = git(["describe", "--exact-match", "--tags", "HEAD"], check=False)
    if r.returncode == 0:
        return r.stdout.strip()
    return None
