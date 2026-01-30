import shutil
import time
from pathlib import Path

import requests
from cookiecutter.main import cookiecutter

from ... import __version__, exceptions, overleaf, paths
from ...logging import get_logger
from ...subproc import get_stdout
from ...zenodo import Zenodo


def setup(slug, cache, overleaf_id, ssh, action_spec, action_version):
    """Set up a new article repo.

    Args:
        slug (str): Repository slug (user/repo).
        cache (bool): If True, enable caching on Zenodo Sandbox.
        overleaf_id (str or NoneType): Overleaf ID of the article.
        ssh (bool): If True, use SSH to clone the repository. Otherwise, use HTTPS.
        action_spec (str or None): Showyourwork version passed to showyourwork-action in
            `.github/workflows/*.yml`
        action_version (str or None):
            Version of the showyourwork-action to use in the workflow.
    """

    # Parse the slug
    user, repo = slug.split("/")
    if Path(repo).exists():
        # Check if it's an empty git repository
        def callback(code, stdout, stderr):
            if code != 0:
                # `git log` errors if there are no commits;
                # that's what we want
                return
            else:
                # There is some history; let's not re-write it
                raise exceptions.ShowyourworkException(
                    f"Directory already exists: {repo}."
                )

        # Require nothing but a `.git` folder and no commit history
        contents = set(list(Path(repo).glob("*"))) - set([Path(repo) / ".git"])
        if len(contents) == 0:
            get_stdout("git log 2> /dev/null", shell=True, cwd=repo, callback=callback)
        else:
            # Refuse to continue
            raise exceptions.ShowyourworkException(f"Directory already exists: {repo}.")

    name = f"@{user}".replace("_", "")

    # Fetch latest showyourwork-action version if not specified
    def get_latest_action_version():
        # Get latest release info
        release_url = "https://api.github.com/repos/showyourwork/showyourwork-action/releases/latest"
        main_url = (
            "https://api.github.com/repos/showyourwork/showyourwork-action/commits/main"
        )
        compare_url = "https://api.github.com/repos/showyourwork/showyourwork-action/compare/{release_tag}...main"
        try:
            release_resp = requests.get(release_url, timeout=5)
            if release_resp.status_code == 200:
                release_data = release_resp.json()
                release_tag = release_data.get("tag_name")
                # release_sha is not needed for comparison logic
                # Get latest commit on main
                main_resp = requests.get(main_url, timeout=5)
                if main_resp.status_code == 200:
                    main_data = main_resp.json()
                    main_sha = main_data.get("sha")
                    # Compare release tag with main
                    if release_tag and main_sha:
                        # Check if main has new commits since release
                        comp_url = compare_url.format(release_tag=release_tag)
                        comp_resp = requests.get(comp_url, timeout=5)
                        if comp_resp.status_code == 200:
                            comp_data = comp_resp.json()
                            ahead_by = comp_data.get("ahead_by", 0)
                            if ahead_by > 0:
                                return main_sha
                        # If no new commits, use release tag
                        return release_tag
        except Exception:
            pass
        # Fallback: just get latest commit SHA from main
        try:
            main_resp = requests.get(main_url, timeout=5)
            if main_resp.status_code == 200:
                main_data = main_resp.json()
                main_sha = main_data.get("sha")
                if main_sha:
                    return main_sha
        except Exception:
            pass
        # Final fallback
        return "v1"

    if not action_version:
        action_version = get_latest_action_version()

    # Set action_spec to match local version if not provided
    if not action_spec:
        DEV_GIT_URL = "git+https://github.com/showyourwork/showyourwork"
        # Get the directory where showyourwork is installed
        showyourwork_dir = str(paths.showyourwork().module)

        # Check git status first to detect uncommitted changes
        has_uncommitted_changes = False
        try:
            git_status = get_stdout(
                "git status --porcelain",
                shell=True,
                cwd=showyourwork_dir,
            ).strip()
            has_uncommitted_changes = bool(git_status)
        except exceptions.CalledProcessError:
            # If git status fails, assume no uncommitted changes
            pass

        # Use same version remotely as locally for consistency
        local_version = __version__

        # If it's a clean release version AND no uncommitted changes, use PyPI
        if not has_uncommitted_changes and not (
            ".dev" in local_version
            or "+" in local_version
            or local_version.endswith(".dirty")
            or "rc" in local_version
        ):
            # Clean release version with clean working tree - use PyPI directly
            action_spec = f"showyourwork=={local_version}"
        else:
            # Development version - need to determine the right git reference
            try:
                # Get current branch
                current_branch = get_stdout(
                    "git branch --show-current",
                    shell=True,
                    cwd=showyourwork_dir,
                ).strip()

                # Use git's tracking branch information to find the right remote
                try:
                    # Get the upstream tracking branch for current branch
                    upstream_branch = get_stdout(
                        "git rev-parse --abbrev-ref --symbolic-full-name @{u}",
                        shell=True,
                        cwd=showyourwork_dir,
                    ).strip()

                    if upstream_branch and "/" in upstream_branch:
                        # Parse remote/branch format (e.g., "upstream/fix-branch")
                        remote_name, remote_branch = upstream_branch.split("/", 1)

                        # Check for uncommitted changes (dirty working tree)
                        try:
                            # Check if working tree is dirty
                            git_status = get_stdout(
                                "git status --porcelain",
                                shell=True,
                                cwd=showyourwork_dir,
                            ).strip()

                            has_uncommitted = bool(git_status)

                            # Check for unpushed commits
                            local_commit = get_stdout(
                                "git rev-parse HEAD",
                                shell=True,
                                cwd=showyourwork_dir,
                            ).strip()

                            remote_commit = get_stdout(
                                f"git rev-parse {upstream_branch}",
                                shell=True,
                                cwd=showyourwork_dir,
                            ).strip()

                            has_unpushed = local_commit != remote_commit

                            # Warn if there are any local changes not reflected
                            # remotely
                            if has_uncommitted or has_unpushed:
                                logger = get_logger()
                                if has_uncommitted and has_unpushed:
                                    logger.warning(
                                        f"Local branch has uncommitted changes AND "
                                        f"unpushed commits. Remote builds will use "
                                        f"{upstream_branch} "
                                        f"(commit {remote_commit[:8]}), "
                                        f"which will differ from your "
                                        f"local working state."
                                    )
                                elif has_uncommitted:
                                    logger.warning(
                                        f"Local branch has uncommitted changes. "
                                        f"Remote builds will use {upstream_branch} "
                                        f"(commit {remote_commit[:8]}), which may "
                                        f"differ from your local working state."
                                    )
                                else:  # has_unpushed
                                    logger.warning(
                                        f"Local branch has unpushed commits. "
                                        f"Remote builds will use {upstream_branch} "
                                        f"(commit {remote_commit[:8]}), which may "
                                        f"differ from your local code."
                                    )
                        except exceptions.CalledProcessError:
                            # If any git commands fail, just proceed
                            pass

                        # Get the URL of this remote
                        remote_url = get_stdout(
                            f"git config --get remote.{remote_name}.url",
                            shell=True,
                            cwd=showyourwork_dir,
                        ).strip()

                        if "github.com" in remote_url:
                            # Parse GitHub URL to get owner/repo
                            if remote_url.startswith("git@"):
                                repo_part = remote_url.split(":")[-1].replace(
                                    ".git", ""
                                )
                            else:
                                repo_part = remote_url.split("github.com/")[-1].replace(
                                    ".git", ""
                                )

                            owner, repo_name = repo_part.split("/")
                            action_spec = (
                                f"git+https://github.com/{owner}/"
                                f"{repo_name}@{remote_branch}"
                            )
                        else:
                            # Not a GitHub remote, fallback
                            action_spec = DEV_GIT_URL
                    else:
                        # No tracking branch set - warn and use main
                        logger = get_logger()
                        logger.warning(
                            f"Branch '{current_branch}' has no remote tracking "
                            f"branch. Using main branch for remote builds. "
                            f"Consider pushing your branch first with: "
                            f"git push -u <remote> {current_branch}"
                        )
                        action_spec = DEV_GIT_URL

                except exceptions.CalledProcessError:
                    # If git tracking fails, fallback to main
                    action_spec = DEV_GIT_URL

            except exceptions.CalledProcessError:
                # If git commands fail, fallback to official main branch
                action_spec = DEV_GIT_URL

    # Create a Zenodo deposit draft for this repo
    if cache:
        deposit_sandbox = Zenodo("sandbox", slug=slug, branch="main")
        cache_sandbox_doi = deposit_sandbox.doi
    else:
        deposit_sandbox = None
        cache_sandbox_doi = ""

    # Set up the repo
    cookiecutter(
        str(paths.showyourwork().cookiecutter),
        no_input=True,
        extra_context={
            "user": user,
            "repo": repo,
            "name": name,
            "showyourwork_version": __version__,
            "cache_sandbox_doi": cache_sandbox_doi,
            "overleaf_id": overleaf_id,
            "year": time.localtime().tm_year,
            "action_spec": action_spec,
            "action_version": action_version,
        },
        overwrite_if_exists=True,
    )

    # Set up git
    try:
        get_stdout("git init -q", shell=True, cwd=repo)
        get_stdout("git add .", shell=True, cwd=repo)
        get_stdout(
            'git commit -q -m "[showyourwork] first commit"',
            shell=True,
            cwd=repo,
        )
        get_stdout("git branch -M main", shell=True, cwd=repo)

        # Set up the remote if it doesn't exist
        def callback(code, stdout, stderr):
            if code != 0:
                if ssh:
                    get_stdout(
                        f"git remote add origin git@github.com:{user}/{repo}.git",
                        shell=True,
                        cwd=repo,
                    )
                else:
                    get_stdout(
                        f"git remote add origin https://github.com/{user}/{repo}.git",
                        shell=True,
                        cwd=repo,
                    )

        get_stdout(
            "git config --get remote.origin.url",
            shell=True,
            cwd=repo,
            callback=callback,
        )

    except Exception as e:
        logger = get_logger()
        if cache:
            logger.error("Deleting cache deposit...")
            deposit_sandbox.delete()
        logger.error(f"Deleting directory {repo}...")
        shutil.rmtree(repo)
        raise e

    # Set up repository on overleaf
    if overleaf_id is not None:
        try:
            overleaf.setup_remote(overleaf_id, path=Path(repo).absolute())
        except Exception as e:
            logger = get_logger()
            if cache:
                logger.error("Deleting cache deposit...")
                deposit_sandbox.delete()
            logger.error(f"Deleting directory {repo}...")
            shutil.rmtree(repo)
            raise e
