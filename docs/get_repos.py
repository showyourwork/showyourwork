import urllib
from urllib.request import Request, urlopen
import json
import os
import sys


def get_repos(
    filename="showyourwork.yml",
    path=".github/workflows",
    maxpages=10,
    exclude_repos=["rodluger"],
):
    """
    Return a list of repos that use showyourwork.

    """
    repos = []
    dates = []
    for page in range(1, maxpages + 1):
        req = Request(
            f"https://api.github.com/search/code?q=path:{path}+filename:{filename}&per_page=100&page={page}"
        )
        req.add_header("Accept", "application/vnd.github.v3+json")
        API_KEY = os.getenv("GH_API_KEY", None)
        if API_KEY is None:
            # We can't authenticate!
            return []
        req.add_header("Authorization", f"token {API_KEY}")
        content = urlopen(req).read()
        content = json.loads(content)
        if len(content["items"]) == 0:
            break
        for res in content["items"]:
            repo = res["repository"]["full_name"]
            if repo not in exclude_repos:
                repos.append(repo)

                # Get the timestamp when it was last pushed to
                req = Request(f"https://api.github.com/repos/{repo}")
                req.add_header("Accept", "application/vnd.github.v3+json")
                req.add_header("Authorization", f"token {API_KEY}")
                content = urlopen(req).read()
                content = json.loads(content)
                date = content["pushed_at"]
                dates.append(date)

    # Sort by the timestamps
    repos = [repo for _, repo in sorted(zip(dates, repos))]
    repos = repos[::-1]

    return repos