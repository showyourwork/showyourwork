import re
import shutil
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Generator, Iterable, NamedTuple, Optional, Tuple

from plumbum import ProcessExecutionError, local

from . import exceptions, logging, paths

OVERLEAF_BLANK_PROJECT_REGEX_TEMPLATE = r"[\n\r\s]+".join(
    [
        r"\\documentclass{article}",
        r"\\usepackage\[utf8\]{inputenc}",
        r"\\title{[^\n{}]+?}",
        r"\\author{[^\n{}]+?}",
        r"\\date{[^\n{}]+?}",
        r"\\begin{document}",
        r"\\maketitle",
        r"\\section{Introduction}",
        r"\\end{document}",
    ]
)

OVERLEAF_BLANK_PROJECT = r"""\documentclass{article}
\usepackage[utf8]{inputenc}

\title{blank project}
\author{Rodrigo Luger}
\date{April 2022}

\begin{document}

\maketitle

\section{Introduction}

\end{document}
"""


class ApplyConflict(exceptions.OverleafError):
    def __init__(self, level: str = "debug"):
        super().__init__(
            "There was an unsolvable conflict when applying a diff. "
            "This generally only occurs when the stored commit hashes "
            "for the local and/or Overleaf repo get out of sync. "
            "To resolve this issue, make sure that the remote Overleaf "
            "documents are all up to date and then run "
            "'showyourwork sync --force' to force acceptance of remote "
            "changes.",
            level=level,
        )


class Conflict(exceptions.OverleafError):
    def __init__(self, message: str, repo: "Repo", level: str = "error"):
        conflicts = repo.git(
            "diff", "--name-only", "--diff-filter=U"
        ).splitlines()
        conflicts = "\n".join(f"- {f}" for f in conflicts)
        super().__init__(
            f"{message}; "
            f"fix the conflicts in '{repo.path}' and mark them as resolved "
            "using 'git add', then run 'showyourwork sync --continue' to finish"
            f"\nconflicts were found in the following files:\n{conflicts}",
            level=level,
        )


class RebaseConflict(Conflict):
    def __init__(self, new_sha: str, repo: "Repo", level: str = "error"):
        self.new_sha = new_sha
        super().__init__(
            "Unable to automatically merge remote changes into local",
            repo,
            level=level,
        )


class MergeConflict(Conflict):
    def __init__(self, repo: "Repo", level: str = "error"):
        super().__init__(
            "Unable to automatically merge local changes into remote",
            repo,
            level=level,
        )


class PushError(exceptions.OverleafError):
    def __init__(self, stdout, stderr, level: str = "error"):
        super().__init__(
            "Overleaf rejected the pushed changes with the following "
            "message:\n\n"
            f"captured stdout:\n{stdout}\n\ncaptured stderr:\n{stderr}",
            level=level,
        )


class MergeOrRebaseResult(NamedTuple):
    changed: bool
    fast_forward: bool
    updated: "Repo"


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

    def git(self, *args: str, with_config: bool = False, **kwargs: str) -> str:
        if with_config:
            args = (
                "-c",
                "core.editor=true",
                "-c",
                "user.name=showyourwork",
                "-c",
                "user.email=showyourwork@showyourwork",
                *args,
            )
        with local.cwd(self.path):
            return local["git"](*args, **kwargs)

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
                    local["git"]("init")

            else:
                local["git"]("clone", "--branch", self.branch, self.url, d)
                if old:
                    with local.cwd(d):
                        local["git"]("reset", "--hard", self.base_sha)

            yield self._replace(path=Path(d))

    def add_and_commit(
        self,
        message: str,
        include: Optional[Iterable[str]] = None,
        exclude: Optional[Iterable[str]] = None,
    ):
        with local.cwd(self.path):
            add_cmd = local["git"]["add"]
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
        self.commit(message)

    def commit(self, message):
        try:
            self.git(
                "commit", "--allow-empty", "-am", message, with_config=True
            )
        except ProcessExecutionError:
            # The first commit may fail if pre-commit is in use
            self.git(
                "commit", "--allow-empty", "-am", message, with_config=True
            )

    def diff(
        self,
        include: Optional[Iterable[str]] = None,
        exclude: Optional[Iterable[str]] = None,
    ):
        with self.checkout_temp(
            old=True
        ) as old_copy, self.checkout_temp() as new_copy:
            with local.cwd(old_copy.path):
                git = local["git"]
                local["rm"]("-rf", ".git", retcode=None)
                git("init", retcode=None)
                old_copy.add_and_commit(
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
            apply_cmd = local["git"]["apply", "--whitespace=fix"]
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
        # sha = self.current_sha()

        # First we try to just apply the diff. If we can, then it's possible to
        # do a fast-forward merge.
        try:
            self._do_apply(diff, exclude=exclude)

        except ApplyConflict:
            if require_fast_forward:
                raise ApplyConflict(level="error")

        else:
            return True, self.base_sha

        # If we get here, we're going to need to merge these divergent histories
        # using a rebase.
        with self.checkout_temp(old=True) as old_copy:
            # Apply the diff starting from the SHA when we last synced
            old_copy._do_apply(diff, exclude=exclude)

            # Make a dummy commit on this old copy
            old_copy.git("checkout", "-b", branch)
            old_copy.add_and_commit(message, exclude=exclude)
            sha = old_copy.current_sha()

            # Rebase the current version onto the patched old one
            with local.cwd(self.path):
                git = local["git"]
                git(
                    "remote",
                    "add",
                    upstream,
                    f"file://{old_copy.path}",
                )
                git("fetch", upstream, branch)
                try:
                    git("rebase", "FETCH_HEAD")

                except ProcessExecutionError:
                    raise RebaseConflict(sha, self, level="debug")

                finally:
                    git("remote", "rm", upstream)
            return False, sha

    def merge_or_rebase(
        self,
        remote_repo: "Repo",
        include: Optional[Iterable[str]] = None,
        exclude: Optional[Iterable[str]] = None,
        require_fast_forward: bool = False,
    ) -> MergeOrRebaseResult:
        """Merge or rebase the changes from another repo into this one"""
        if self.is_dirty():
            raise exceptions.OverleafError(
                "Local project is dirty; please commit or stash before syncing"
            )
        diff = remote_repo.diff(include=include, exclude=exclude)
        if diff.strip() == "":
            return MergeOrRebaseResult(
                changed=False, fast_forward=True, updated=self
            )
        ff, sha = self.apply(
            diff, exclude=exclude, require_fast_forward=require_fast_forward
        )
        return MergeOrRebaseResult(
            changed=True, fast_forward=ff, updated=self._replace(base_sha=sha)
        )


class Overleaf(NamedTuple):
    local: Repo
    remote: Repo

    @classmethod
    def from_project_id(
        cls,
        project_id: str,
        path: Optional[str] = None,
        local_sha: Optional[str] = None,
        remote_sha: Optional[str] = None,
    ) -> "Overleaf":
        return cls.from_url(
            f"https://git.overleaf.com/{project_id}",
            path=path,
            local_sha=local_sha,
            remote_sha=remote_sha,
        )

    @classmethod
    def from_url(
        cls,
        url: str,
        path: Optional[str] = None,
        local_sha: Optional[str] = None,
        remote_sha: Optional[str] = None,
    ) -> "Overleaf":
        user_paths = paths.user(path=path)
        with local.cwd(user_paths.repo):
            local_branch = local["git"](
                "rev-parse", "--abbrev-ref", "HEAD"
            ).strip()
        return cls(
            local=Repo(
                url=f"file://{user_paths.repo}",
                branch=local_branch,
                path=user_paths.repo,
                base_sha=local_sha,
                subdirectory=str(user_paths.tex.relative_to(user_paths.repo)),
            ),
            remote=Repo(
                url=url,
                branch="master",
                path=user_paths.overleaf,
                base_sha=remote_sha,
            ),
        )

    def to_current(self) -> "Overleaf":
        return self._replace(
            local=self.local._replace(base_sha=self.local.current_sha()),
            remote=self.remote._replace(base_sha=self.remote.current_sha()),
        )

    def init_remote(self):
        if self.remote.path.exists():
            shutil.rmtree(self.remote.path)
        self.remote.path.mkdir(exist_ok=True)
        self.remote.git("init")
        self.remote.git("branch", "-M", self.remote.branch)
        self.remote.git("pull", self.remote.url)

    def push_remote(self):
        self.remote.git("push", self.remote.url, self.remote.branch)

    def setup_remote(self, force: bool = False):
        self.init_remote()

        # If we're not forcing the overwrite, check to make sure that the
        # project is blank.
        if not force:
            files = [f for f in self.remote.path.glob("*") if f.name != ".git"]
            if len(files) != 1 or not re.match(
                OVERLEAF_BLANK_PROJECT_REGEX_TEMPLATE, open(files[0]).read()
            ):
                raise exceptions.OverleafError(
                    "Overleaf repository is not a blank project. "
                    "Use force=True to overwrite."
                )

        # Remove all the files in the remote
        self.remote.git("rm", "-r", "*")

        # Copy over the local files
        files = self.local.git(
            "ls-tree",
            "-r",
            "--name-only",
            "HEAD",
            self.local.subdirectory,
        ).splitlines()
        for f in files:
            f = Path(f)
            base_file = f.relative_to(self.local.subdirectory)
            local = self.local.path / f
            remote = self.remote.path / base_file
            if local.is_dir():
                if remote.exists():
                    shutil.rmtree(remote)
                shutil.copytree(local, remote)
            else:
                shutil.copy(local, remote)
            self.remote.git("add", str(base_file))

        # Commit the changes
        self.remote.add_and_commit("[showyourwork] Initial setup")
        self.push_remote()

        # Save the new base SHA
        return self.to_current()

    @property
    def logger(self):
        return logging.get_logger()

    @property
    def _sync_stage_path(self) -> Path:
        return paths.user(path=self.local.path).sync / "stage"

    def _save_sync_stage(self, stage: str):
        with open(self._sync_stage_path, "w") as f:
            f.write(stage)

    def _load_sync_stage(self) -> str:
        if not self._sync_stage_path.exists():
            raise exceptions.OverleafError("No sync is currently in progress")
        with open(self._sync_stage_path, "r") as f:
            return f.read()

    @property
    def _rebase_lock_path(self) -> Path:
        return paths.user(path=self.local.path).sync / "rebase.lock"

    def _save_rebase_lock(self, sha: str):
        with open(self._rebase_lock_path, "w") as f:
            f.write(sha)

    def _load_rebase_lock(self):
        if not self._rebase_lock_path.exists():
            raise exceptions.OverleafError(
                "No rebase is currently in progress"
            )
        with open(self._rebase_lock_path, "r") as f:
            return f.read()

    def sync(self) -> "Overleaf":
        if self._sync_stage_path.exists():
            stage = self._load_sync_stage()
            if stage == "from remote":
                raise exceptions.OverleafError(
                    "A sync is already in progress; you may need to fix rebase "
                    "conflicts in your local directory and then continue the sync"
                )
            elif stage == "to remote":
                raise exceptions.OverleafError(
                    "A sync is already in progress; you may need to fix merge "
                    f"conflicts in your {paths.user(path=self.local.path).overleaf} "
                    "and then continue the sync"
                )
            raise exceptions.OverleafError("Sync already in progress")

        self.init_remote()

        self._save_sync_stage("from remote")
        new_self = self.sync_from_remote()

        new_self._save_sync_stage("to remote")
        new_self = new_self.sync_to_remote()

        new_self._save_sync_stage.unlink(missing_ok=True)

        return new_self

    def continue_sync(self) -> "Overleaf":
        if not self._sync_stage_path.exists():
            raise exceptions.OverleafError("No sync is currently in progress")

        stage = self._load_sync_stage()
        if stage == "from remote":
            new_self = self.continue_sync_from_remote()
            new_self._save_sync_stage("to remote")
            new_self = new_self.sync_to_remote()
        elif stage == "to remote":
            new_self = self.continue_sync_to_remote()
        else:
            raise exceptions.OverleafError(
                f"Unrecognized sync stage name: {stage}"
            )

        new_self._save_sync_stage.unlink(missing_ok=True)
        return new_self

    def sync_from_remote(self) -> "Overleaf":
        self.logger.info("Pulling remote changes from Overleaf")
        try:
            result = self.local.merge_or_rebase(self.remote)
        except RebaseConflict as e:
            new_sha = e.new_sha
            self._save_rebase_lock(new_sha)
            raise RebaseConflict(new_sha, self.local, level="error")
        if not result.changed:
            self.logger.info("No remote changes found")
        elif result.fast_forward:
            self.logger.info(
                "Automatically merged remote changes; "
                "don't forget to push your local changes"
            )
        else:
            self.logger.info(
                "Automatically rebased against remote changes; "
                "you'll need to force push your local changes: "
                f"git push --force origin {self.local.branch}"
            )
        return self.finish_sync_from_remote(result.updated.base_sha)

    def abort_sync_from_remote(self):
        self._load_rebase_lock()
        self.logger.info("Aborting rebase against remote changes")
        self.local.git("rebase", "--abort")
        self._rebase_lock_path.unlink(missing_ok=True)

    def continue_sync_from_remote(self) -> "Overleaf":
        new_sha = self._load_rebase_lock()
        self.logger.info("Continuing rebase against remote changes")
        try:
            self.local.git("rebase", "--continue", with_config=True)
        except ProcessExecutionError as e:
            raise RebaseConflict(new_sha, self.local)
        self.logger.info(
            "Finished rebasing against remote changes; "
            "you'll need to force push your local changes: "
            f"git push --force origin {self.local.branch}"
        )
        return self.finish_sync_from_remote(new_sha)

    def finish_sync_from_remote(self, new_sha: str) -> "Overleaf":
        self._rebase_lock_path.unlink(missing_ok=True)
        return self._replace(local=self.local._replace(base_sha=new_sha))

    def sync_to_remote(self) -> "Overleaf":
        self.logger.info("Pushing local changes to Overleaf")
        self.remote.merge_or_rebase(self.local, require_fast_forward=True)
        return self.finish_sync_to_remote()

    def continue_sync_to_remote(self) -> "Overleaf":
        self.logger.info("Continuing to push local changes to Overleaf")
        return self.finish_sync_to_remote()

    def finish_sync_to_remote(self) -> "Overleaf":
        self.remote.add_and_commit("[showyourwork] Updating overleaf")
        try:
            self.push_remote()
        except ProcessExecutionError:
            # Try to pull and automatically merge; if this fails, the user will
            # need to fix the conflicts.
            try:
                self.remote.git(
                    "pull",
                    "--no-rebase",
                    self.remote.url,
                    self.remote.branch,
                    with_config=True,
                )
            except ProcessExecutionError as e:
                raise MergeConflict(self.remote)

            # If this still fails, something else went wrong and the user will
            # need to act on the information provided
            try:
                self.push_remote()
            except ProcessExecutionError as e:
                raise PushError(e.stdout, e.stderr)

        self.logger.info("Finished pushing local changes to Overleaf")
        return self.to_current()
