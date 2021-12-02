import urllib
from urllib.request import Request, urlopen
import json
import os
import sys
import re


def get_commit_count(repo, API_KEY):
    """
    Return the number of commits to a project.

    https://stackoverflow.com/a/55873469

    """

    req = Request(f"https://api.github.com/repos/{repo}/commits?per_page=1")
    req.add_header("Accept", "application/vnd.github.v3+json")
    req.add_header("Authorization", f"token {API_KEY}")
    resp = urlopen(req)
    content = resp.read()
    content = json.loads(content)
    if len(content):
        commits = int(
            resp.headers["Link"]
            .split(",")[1]
            .split("per_page=1&page=")[1]
            .split(">")[0]
        )
    else:
        commits = 0
    return commits


def get_repos(
    filename="showyourwork.yml",
    path=".github/workflows",
    maxpages=10,
    exclude_repos=[
        "rodluger/showyourwork-template",
        "rodluger/showyourwork-sandbox",
        "rodluger/showyourwork-example-dev",
        "gusbeane/fdbk_eos-temp",
    ],
):
    """
    Return a list of repos that use showyourwork.

    """
    # Get curated list of projects that use showyourwork
    with open("projects.json", "r") as f:
        projects = json.load(f)
    for project in projects:
        projects[project]["date"] = projects[project].get("date", "")

    # Get all the showyourwork tags (versions)
    API_KEY = os.getenv("GH_API_KEY", None)
    if API_KEY is None:
        print("ERROR: Can't authenticate git.")
        return []
    req = Request(f"https://api.github.com/repos/rodluger/showyourwork/tags")
    req.add_header("Accept", "application/vnd.github.v3+json")
    req.add_header("Authorization", f"token {API_KEY}")
    content = urlopen(req).read()
    content = json.loads(content)
    versions = {version["commit"]["sha"]: version["name"] for version in content}

    # Now do an API search for all repos that use showyourwork
    for page in range(1, maxpages + 1):
        req = Request(
            f"https://api.github.com/search/code?q=path:{path}+filename:{filename}&per_page=100&page={page}"
        )
        req.add_header("Accept", "application/vnd.github.v3+json")
        req.add_header("Authorization", f"token {API_KEY}")
        try:
            content = urlopen(req).read()
        except urllib.error.HTTPError as e:
            break
        content = json.loads(content)
        if len(content["items"]) == 0:
            break
        for res in content["items"]:
            repo = res["repository"]["full_name"]
            if repo not in exclude_repos:

                # Get the timestamp when it was last pushed to
                req = Request(f"https://api.github.com/repos/{repo}")
                req.add_header("Accept", "application/vnd.github.v3+json")
                req.add_header("Authorization", f"token {API_KEY}")
                content = urlopen(req).read()
                content = json.loads(content)
                date = content["pushed_at"]

                # Get the SHA or version for the showyourwork submodule
                req = Request(f"https://api.github.com/repos/{repo}/git/trees/main")
                req.add_header("Accept", "application/vnd.github.v3+json")
                req.add_header("Authorization", f"token {API_KEY}")
                content = urlopen(req).read()
                content = json.loads(content)
                try:
                    sha = [c for c in content["tree"] if c["path"] == "showyourwork"][
                        0
                    ]["sha"]
                    version = versions.get(sha, sha[:7])
                except:
                    version = "unknown"

                # Get the number of commits
                try:
                    commits = get_commit_count(repo, API_KEY)
                except:
                    commits = "N/A"

                # Assemble the metadata
                if repo in projects:
                    projects[repo]["version"] = version
                    projects[repo]["commits"] = commits
                    projects[repo]["date"] = date
                else:
                    projects[repo] = {
                        "title": "",
                        "authors": [],
                        "field": "Uncategorized",
                        "version": version,
                        "commits": commits,
                        "date": date,
                    }

    return projects