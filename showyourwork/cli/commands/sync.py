
from ... import logging, paths
from ...config import edit_yaml
from ...git import get_repo_branch
from ...overleaf2 import Overleaf


def overleaf_sync(init=False, force=False):
    branch = get_repo_branch()
    with edit_yaml("overleaf.yml") as config:
        version = config["overleaf"].get("version", 1)
        assert version == 2

        url = config["overleaf"].get(
            "url", f"https://git.overleaf.com/{config['overleaf']['id']}"
        )
        local_sha = config["overleaf"].get("local_sha", None)
        remote_sha = config["overleaf"].get("remote_sha", None)
