from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest
from plumbum import local

from showyourwork import exceptions
from showyourwork.overleaf2 import (
    OVERLEAF_BLANK_PROJECT,
    MergeConflict,
    Overleaf,
    RebaseConflict,
    Repo,
)

git = local["git"]


def _new_repo(path, branch):
    with local.cwd(path):
        git("init", "-b", branch)
        git("add", ".")
        git("config", "user.name", "test")
        git("config", "user.email", "test")
        git("commit", "--allow-empty", "-am", "initial", retcode=None)


@contextmanager
def overleaf_like(empty=False, blank=False):
    with TemporaryDirectory() as d:
        if blank:
            p = Path(d)
            with open(p / "main.tex", "w") as f:
                f.write(OVERLEAF_BLANK_PROJECT)
        elif not empty:
            p = Path(d)
            with open(p / "ms.tex", "w") as f:
                f.write(OVERLEAF_BLANK_PROJECT)
            with open(p / "bib.bib", "w") as f:
                f.write("This is bib.bib\n")
            (p / "figures").mkdir(parents=True)
            with open(p / "figures" / "figure.png", "w") as f:
                pass

        _new_repo(d, "master")
        repo = Repo(url=f"file://{d}", branch="master", path=Path(d))
        repo.git("config", "receive.denyCurrentBranch", "ignore")
        yield repo._replace(base_sha=repo.current_sha())


@contextmanager
def syw_like(empty=False):
    with TemporaryDirectory() as d:
        if not empty:
            p = Path(d) / "src" / "tex"
            p.mkdir(parents=True)
            with open(Path(d) / ".gitignore", "w") as f:
                f.write(".showyourwork\n")
            with open(p / "ms.tex", "w") as f:
                f.write(OVERLEAF_BLANK_PROJECT)
            with open(p / "bib.bib", "w") as f:
                f.write("This is bib.bib\n")

        _new_repo(d, "main")
        repo = Repo(
            url=f"file://{d}",
            branch="main",
            path=Path(d),
            subdirectory="src/tex",
        )
        yield repo._replace(base_sha=repo.current_sha())


@contextmanager
def repo_pair(order, syw_empty=False, ovl_empty=False, blank=False):
    with syw_like(empty=syw_empty) as local_repo, overleaf_like(
        empty=ovl_empty, blank=blank
    ) as remote_repo:
        if order:
            yield local_repo, remote_repo
        else:
            yield remote_repo, local_repo


@pytest.mark.parametrize("order", [True, False])
def test_dirty_repo(order):
    with repo_pair(order=order) as (syw, ovl):
        with open(syw.source_path / "ms.tex", "a") as f:
            f.write("This is an additional line\n")
        with pytest.raises(exceptions.OverleafError):
            syw.merge_or_rebase(ovl)


@pytest.mark.parametrize("order", [True, False])
def test_fast_forward(order):
    with repo_pair(order) as (syw, ovl):
        fn = ovl.source_path / "ms.tex"
        open(fn, "a").write("This is an additional line\n")
        expected = open(fn, "r").read()
        ovl.git("commit", "-am", "additional line")
        syw.merge_or_rebase(ovl)
        assert open(syw.source_path / "ms.tex", "r").read() == expected


@pytest.mark.parametrize("order", [True, False])
@pytest.mark.parametrize("exclude", ["figures/*", "figures/*.png"])
def test_exclude(order, exclude):
    with repo_pair(order, syw_empty=order, ovl_empty=not order) as (syw, ovl):
        syw.merge_or_rebase(ovl._replace(base_sha=None), exclude=[exclude])
        assert (
            open(syw.source_path / "ms.tex", "r").read()
            == open(ovl.source_path / "ms.tex", "r").read()
        )
        assert not (syw.source_path / "figures").exists()


@pytest.mark.parametrize("order", [True, False])
def test_merge_same_file(order):
    with repo_pair(order) as (syw, ovl):
        ovl_fn = ovl.source_path / "ms.tex"
        open(ovl_fn, "a").write("This is an additional line\n")
        ovl.git("commit", "-am", "additional line: ovl")

        syw_fn = syw.source_path / "ms.tex"
        current = open(syw_fn, "r").read()
        open(syw_fn, "w").write("This is a new first line\n" + current)
        syw.git("commit", "-am", "additional line: syw")

        expected = open(syw_fn, "r").read() + open(ovl_fn, "r").readlines()[-1]

        ff, _ = syw.merge_or_rebase(ovl)
        assert ff
        assert open(syw_fn, "r").read() == expected


@pytest.mark.parametrize("order", [True, False])
def test_merge_different_files(order):
    with repo_pair(order) as (syw, ovl):
        ovl_fn = ovl.source_path / "ms.tex"
        open(ovl_fn, "a").write("This is an additional line\n")
        ovl.git("commit", "-am", "additional line: ovl")
        expected_ms = open(ovl_fn, "r").read()

        syw_fn = syw.source_path / "bib.bib"
        open(syw_fn, "a").write("This is an additional line in the bib\n")
        syw.git("commit", "-am", "additional line: syw")
        expected_bib = open(syw_fn, "r").read()

        ff, _ = syw.merge_or_rebase(ovl)
        assert ff
        assert open(syw.source_path / "ms.tex", "r").read() == expected_ms
        assert open(syw.source_path / "bib.bib", "r").read() == expected_bib


@pytest.mark.parametrize("order", [True, False])
def test_conflict(order):
    with repo_pair(order) as (syw, ovl):
        ovl_fn = ovl.source_path / "ms.tex"
        open(ovl_fn, "a").write("This is an additional line\n")
        ovl.git("commit", "-am", "additional line: ovl")

        syw_fn = syw.source_path / "ms.tex"
        open(syw_fn, "a").write("This is a different additional line\n")
        syw.git("commit", "-am", "additional line: syw")

        with pytest.raises(RebaseConflict):
            syw.merge_or_rebase(ovl)

        expected = OVERLEAF_BLANK_PROJECT + "This is an additional line\n"
        open(syw_fn, "w").write(expected)
        syw.git("add", syw_fn)
        syw.git("rebase", "--continue")
        assert open(syw_fn, "r").read() == expected


@contextmanager
def setup_overleaf(**kwargs):
    kwargs["blank"] = kwargs.get("blank", True)
    with repo_pair(True, **kwargs) as (syw, ovl):
        overleaf = Overleaf.from_url(ovl.url, path=syw.path)
        overleaf = overleaf.setup_remote()
        ovl.git("reset", "--hard", "HEAD")
        yield overleaf, ovl


def sync(overleaf, source_repo, reinit=True):
    # This re-pulls the remote, incorporating the changes we just made. If
    # we didn't do this, it would be as if these changes were made during
    # the sync process
    if reinit:
        overleaf.init_remote()

    overleaf = overleaf.sync_from_remote()
    overleaf = overleaf.sync_to_remote()

    # Sync the pushed changes
    source_repo.git("reset", "--hard", "HEAD")

    return overleaf


def test_setup():
    with setup_overleaf() as (overleaf, _):
        assert not (overleaf.remote.source_path / "main.tex").exists()
        assert (
            open(overleaf.remote.source_path / "ms.tex").read()
            == open(overleaf.local.source_path / "ms.tex").read()
        )


def test_setup_not_blank():
    with pytest.raises(exceptions.OverleafError):
        with setup_overleaf(blank=False):
            pass


def test_sync_remote_changes():
    with setup_overleaf() as (overleaf, source_repo):
        with open(source_repo.source_path / "ms.tex", "a") as f:
            f.write("This is an additional line\n")
        source_repo.add_and_commit("Adding an additional line")
        overleaf = sync(overleaf, source_repo)
        expect = open(overleaf.local.source_path / "ms.tex").read()
        assert open(source_repo.source_path / "ms.tex").read() == expect
        assert expect.endswith("This is an additional line\n")


def test_sync_local_changes():
    with setup_overleaf() as (overleaf, source_repo):
        with open(overleaf.local.source_path / "ms.tex", "a") as f:
            f.write("This is an additional line\n")
        overleaf.local.add_and_commit("Adding an additional line")
        overleaf = sync(overleaf, source_repo)
        expect = open(overleaf.local.source_path / "ms.tex").read()
        assert open(source_repo.source_path / "ms.tex").read() == expect
        assert expect.endswith("This is an additional line\n")


def test_sync_both_changes():
    with setup_overleaf() as (overleaf, source_repo):
        with open(source_repo.source_path / "ms.tex", "a") as f:
            f.write("This is an additional line\n")
        source_repo.add_and_commit("Adding an additional line")

        data = open(overleaf.local.source_path / "ms.tex").read()
        with open(overleaf.local.source_path / "ms.tex", "w") as f:
            f.write("This is a new first line\n" + data)
        overleaf.local.add_and_commit("Adding a different additional line")

        overleaf = sync(overleaf, source_repo)
        expect = open(overleaf.local.source_path / "ms.tex").read()
        assert open(source_repo.source_path / "ms.tex").read() == expect
        assert expect.startswith("This is a new first line\n")
        assert expect.endswith("This is an additional line\n")


def test_sync_rebase_abort():
    with setup_overleaf() as (overleaf, source_repo):
        with open(source_repo.source_path / "ms.tex", "a") as f:
            f.write("This is an additional line\n")
        source_repo.add_and_commit("Adding an additional line")

        with open(overleaf.local.source_path / "ms.tex", "a") as f:
            f.write("This is a different additional line\n")
        overleaf.local.add_and_commit("Adding a different additional line")
        expect = open(overleaf.local.source_path / "ms.tex").read()

        with pytest.raises(RebaseConflict):
            sync(overleaf, source_repo)

        overleaf.abort_sync_from_remote()
        assert open(overleaf.local.source_path / "ms.tex").read() == expect


def test_sync_rebase_merge():
    with setup_overleaf() as (overleaf, source_repo):
        with open(source_repo.source_path / "ms.tex", "a") as f:
            f.write("This is an additional line\n")
        source_repo.add_and_commit("Adding an additional line")

        baseline = open(overleaf.local.source_path / "ms.tex").read()
        with open(overleaf.local.source_path / "ms.tex", "a") as f:
            f.write("This is a different additional line\n")
        overleaf.local.add_and_commit("Adding a different additional line")

        with pytest.raises(RebaseConflict):
            sync(overleaf, source_repo)

        # Fix the conflict
        with open(overleaf.local.source_path / "ms.tex", "w") as f:
            f.write(baseline + "This is a rebased additional line\n")
        overleaf.local.git("add", overleaf.local.source_path / "ms.tex")

        overleaf = overleaf.continue_sync_from_remote()
        overleaf = overleaf.sync_to_remote()
        source_repo.git("reset", "--hard", "HEAD")

        expect = open(overleaf.local.source_path / "ms.tex").read()
        assert open(source_repo.source_path / "ms.tex").read() == expect
        assert expect.endswith("This is a rebased additional line\n")


def test_remote_changes_during_sync():
    with setup_overleaf() as (overleaf, source_repo):
        initial = open(source_repo.source_path / "ms.tex", "r").read()
        with open(source_repo.source_path / "ms.tex", "a") as f:
            f.write("This is an additional line\n")
        source_repo.add_and_commit("Adding an additional line")
        overleaf.init_remote()

        # Add a change that happens during the sync
        expect = initial + "This is a different additional line\n"
        with open(source_repo.source_path / "ms.tex", "w") as f:
            f.write(expect)
        source_repo.add_and_commit("Adding a different line")

        overleaf = sync(overleaf, source_repo, reinit=False)
        assert open(source_repo.source_path / "ms.tex").read() == expect


def test_remote_changes_during_sync_with_conflict():
    with setup_overleaf() as (overleaf, source_repo):
        initial = open(overleaf.local.source_path / "ms.tex", "r").read()
        with open(overleaf.local.source_path / "ms.tex", "a") as f:
            f.write("This is an additional line\n")
        overleaf.local.add_and_commit("Adding an additional line")
        overleaf = overleaf.sync_from_remote()

        # Add a change that happens during the sync
        expect = initial + "This is a different additional line\n"
        with open(source_repo.source_path / "ms.tex", "w") as f:
            f.write(expect)
        source_repo.add_and_commit("Adding a different line")

        # This should fail with a merge conflict
        with pytest.raises(MergeConflict):
            overleaf.sync_to_remote()

        # Fix the merge conflict and finish the sync
        with open(overleaf.local.source_path / "ms.tex", "w") as f:
            f.write(expect)
        overleaf.local.git("add", overleaf.local.source_path / "ms.tex")
        overleaf = overleaf.continue_sync_to_remote()

        assert open(source_repo.source_path / "ms.tex").read() == expect
