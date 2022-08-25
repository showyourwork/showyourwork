"""
Implements the logic for merging two disparate git histories. This is used for
Overleaf support, but implemented separately to make testing and maintenance
easier.
"""

from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Generator, Iterable, NamedTuple, Optional, Tuple

from plumbum import ProcessExecutionError, local

from . import exceptions

git = local["git"]
rm = local["rm"]


class MergeError(exceptions.ShowyourworkException):
    pass


class ApplyConflict(exceptions.ShowyourworkException):
    def __init__(self):
        super().__init__(
            "Cannot apply diff. This generally means that you made local "
            "(or remote) changes to a file that is owned by remote (or "
            "local); re-run `sync` with the argument `--force` to "
            "overwrite changes"
        )


class RebaseConflict(exceptions.ShowyourworkException):
    def __init__(self, stdout, stderr):
        super().__init__(
            "Unable to automatically merge histories; "
            "fix conflicts and continue the rebase as described below:\n\n"
            f"{stdout}\n{stderr}"
        )


class Repo(NamedTuple):
    url: str
    branch: str
    path: Path
    base_sha: Optional[str] = None
    subdirectory: Optional[str] = None

    @property
    def source_path(self) -> Path:
        if self.subdirectory is None:
            return self.path
        else:
            return self.path / self.subdirectory

    def git(self, *args: str, **kwargs: str) -> str:
        with local.cwd(self.path):
            return git(*args, **kwargs)

    def current_sha(self) -> str:
        return self.git("rev-parse", "HEAD").strip()

    def is_dirty(self) -> bool:
        return bool(self.git("status", "--porcelain").strip())

    @contextmanager
    def checkout_temp(
        self, old: bool = False
    ) -> Generator["Repo", None, None]:
        with TemporaryDirectory() as d:
            # If there's no base commit, we want the whole history
            if old and self.base_sha is None:
                with local.cwd(d):
                    git("init")

            else:
                git("clone", "--branch", self.branch, self.url, d)
                if old:
                    with local.cwd(d):
                        git("reset", "--hard", self.base_sha)

            yield self._replace(path=Path(d))

    def _add_and_commit(
        self,
        message: str,
        include: Optional[Iterable[str]] = None,
        exclude: Optional[Iterable[str]] = None,
    ):
        with local.cwd(self.path):
            add_cmd = git["add"]
            if include is not None:
                assert not exclude
                if self.subdirectory is None:
                    for p in include:
                        add_cmd = add_cmd[f"{p}"]
                else:
                    for p in include:
                        add_cmd = add_cmd[f"{Path(self.subdirectory) / p}"]
            elif exclude is None:
                add_cmd = add_cmd["."]
            else:
                add_cmd = add_cmd["--", "."]
                if self.subdirectory is None:
                    for p in exclude:
                        add_cmd = add_cmd[f":!{p}"]
                else:
                    for p in exclude:
                        add_cmd = add_cmd[f":!{Path(self.subdirectory) / p}"]
            add_cmd()
            git("config", "user.name", "showyourwork")
            git("config", "user.email", "showyourwork@showyourwork")
            try:
                git("commit", "--allow-empty", "-am", message)
            except ProcessExecutionError:
                # The first commit may fail if pre-commit is in use
                git("commit", "--allow-empty", "-am", message)
            git("config", "--unset", "user.name")
            git("config", "--unset", "user.email")

    def diff(
        self,
        include: Optional[Iterable[str]] = None,
        exclude: Optional[Iterable[str]] = None,
    ):
        with self.checkout_temp(
            old=True
        ) as old_copy, self.checkout_temp() as new_copy:
            with local.cwd(old_copy.path):
                rm("-rf", ".git", retcode=None)
                git("init", retcode=None)
                old_copy._add_and_commit(
                    "dummy commit", include=include, exclude=exclude
                )
                git(
                    "remote",
                    "add",
                    "real_dst",
                    f"file://{new_copy.path}",
                )
                git("fetch", "--depth=1", "real_dst", "HEAD")
                if self.subdirectory is None:
                    return git(
                        "diff-tree",
                        "--unified=1",
                        "HEAD...FETCH_HEAD",
                        "--inter-hunk-context=-1",
                    )
                else:
                    return git(
                        "diff-tree",
                        "--unified=1",
                        "HEAD...FETCH_HEAD",
                        "--inter-hunk-context=-1",
                        f"--relative={self.subdirectory}",
                    )

    def _do_apply(
        self, diff: str, exclude: Optional[Iterable[str]] = None
    ) -> None:
        with local.cwd(self.path):
            apply_cmd = git["apply", "--reject", "--whitespace=fix"]
            if self.subdirectory is not None:
                (self.path / self.subdirectory).mkdir(
                    parents=True, exist_ok=True
                )
                apply_cmd = apply_cmd[f"--directory={self.subdirectory}"]
            if exclude:
                if self.subdirectory is None:
                    for p in exclude:
                        apply_cmd = apply_cmd[f"--exclude={p}"]
                else:
                    for p in exclude:
                        apply_cmd = apply_cmd[
                            f"--exclude={Path(self.subdirectory) / p}"
                        ]
            try:
                (apply_cmd << diff)()
            except ProcessExecutionError:
                raise ApplyConflict()

    def apply(
        self,
        diff: str,
        message: str = "applying diff",
        upstream: str = "showyourwork_auto_upstream",
        branch: str = "showyourwork_auto_branch",
        exclude: Optional[Iterable[str]] = None,
        require_fast_forward: bool = False,
    ) -> Tuple[bool, str]:
        sha = self.current_sha()

        # First we try to just apply the diff. If we can, then it's possible to
        # do a fast-forward merge.
        try:
            self._do_apply(diff, exclude=exclude)

        except ApplyConflict:
            if require_fast_forward:
                raise ApplyConflict()

        else:
            return True, sha

        # If we get here, we're going to need to merge these divergent histories
        # using a rebase.
        with self.checkout_temp(old=True) as old_copy:
            # Apply the diff starting from the SHA when we last synced
            old_copy._do_apply(diff, exclude=exclude)

            # Make a dummy commit on this old copy
            old_copy.git("checkout", "-b", branch)
            old_copy._add_and_commit(message, exclude=exclude)
            sha = old_copy.current_sha()

            # Rebase the current version onto the patched old one
            with local.cwd(self.path):
                git(
                    "remote",
                    "add",
                    upstream,
                    f"file://{old_copy.path}",
                )
                git("fetch", upstream, branch)
                try:
                    git("rebase", "FETCH_HEAD")

                except ProcessExecutionError as e:
                    raise RebaseConflict(e.stdout, e.stderr)

                finally:
                    git("remote", "rm", upstream)
            return False, sha

    def merge_or_rebase(
        self,
        remote_repo: "Repo",
        include: Optional[Iterable[str]] = None,
        exclude: Optional[Iterable[str]] = None,
        require_fast_forward: bool = False,
    ) -> Tuple[bool, "Repo"]:
        """Merge or rebase the changes from another repo into this one

        If possible this will
        """
        if self.is_dirty():
            raise exceptions.OverleafError(
                "Local project is dirty; please commit or stash before syncing"
            )
        diff = remote_repo.diff(include=include, exclude=exclude)
        ff, sha = self.apply(
            diff, exclude=exclude, require_fast_forward=require_fast_forward
        )
        return ff, self._replace(base_sha=sha)
