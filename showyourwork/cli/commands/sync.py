from ... import logging, paths
from ...config import edit_yaml
from ...overleaf2 import Overleaf


def overleaf_sync(continue_=False, init=False, force=False):
    with edit_yaml("overleaf.yml") as config:
        version = config["overleaf"].get("version", 1)
        assert version == 2

        url = config["overleaf"].get(
            "url", f"https://git.overleaf.com/{config['overleaf']['id']}"
        )
        local_sha = config["overleaf"].get("local_sha", None)
        remote_sha = config["overleaf"].get("remote_sha", None)

    overleaf = Overleaf.from_url(
        url, local_sha=local_sha, remote_sha=remote_sha
    )

    if init:
        overleaf.setup_remote(force=force)
        return

    if continue_:
        overleaf = overleaf.continue_sync()

    else:
        overleaf = overleaf.sync(force=force)
