import urllib
from urllib.request import Request, urlopen
import json
import os
import sys
import re
import time


def get_commit_count(repo, API_KEY):
    """
    Return the number of commits to a project.

    https://stackoverflow.com/a/55873469

    """

    req = Request(f"https://api.github.com/repos/{repo}/commits?per_page=1")
    req.add_header("Accept", "application/vnd.github.v3+json")
    req.add_header("Authorization", f"token {API_KEY}")
    resp = urlopen(req)
    try:
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
    except:
        commits = "N/A"
    return commits


def get_date(repo, API_KEY):
    """
    Get the timestamp when a repo was last pushed to.

    """
    req = Request(f"https://api.github.com/repos/{repo}")
    req.add_header("Accept", "application/vnd.github.v3+json")
    req.add_header("Authorization", f"token {API_KEY}")
    req.add_header("User-Agent", "request")
    try:
        content = urlopen(req).read()
        content = json.loads(content)
        date = content["pushed_at"]
    except:
        date = "??"
    return date


def get_version(repo, versions, API_KEY):
    """
    Get the SHA or version for the showyourwork submodule in the repo.

    """
    req = Request(f"https://api.github.com/repos/{repo}/git/trees/main")
    req.add_header("Accept", "application/vnd.github.v3+json")
    req.add_header("Authorization", f"token {API_KEY}")
    req.add_header("User-Agent", "request")
    try:
        content = urlopen(req).read()
        content = json.loads(content)
        sha = [c for c in content["tree"] if c["path"] == "showyourwork"][0]["sha"]
        version = versions.get(sha, sha[:7])
    except:
        version = "unknown"
    return version


def get_repos(
    filename="showyourwork.yml",
    path=".github/workflows",
    maxpages=10,
    maxtries=3,
    sleep_time=5.0,
    exclude_repos=[
        "rodluger/showyourwork-template",
        "rodluger/showyourwork-sandbox",
        "rodluger/showyourwork-example-dev",
        "gusbeane/fdbk_eos-temp",
        "LBJ-Wade/showyourwork-template",
    ],
):
    """
    Return a list of repos that use showyourwork.

    """
    # Get all the showyourwork tags (versions)
    API_KEY = os.getenv("GH_API_KEY", None)
    if API_KEY is None:
        print("ERROR: Can't authenticate git.")
        return []
    req = Request(f"https://api.github.com/repos/rodluger/showyourwork/tags")
    req.add_header("Accept", "application/vnd.github.v3+json")
    req.add_header("Authorization", f"token {API_KEY}")
    req.add_header("User-Agent", "request")
    content = urlopen(req).read()
    content = json.loads(content)
    versions = {version["commit"]["sha"]: version["name"] for version in content}

    # Get curated list of projects that use showyourwork
    with open("projects.json", "r") as f:
        projects = json.load(f)
    for project in projects:
        projects[project]["date"] = projects[project].get("date", "")
        projects[project]["doi"] = projects[project].get("doi", "N/A")
        projects[project]["url"] = projects[project].get("url", "")
        projects[project]["version"] = get_version(project, versions, API_KEY)
        projects[project]["commits"] = get_commit_count(project, API_KEY)
        projects[project]["date"] = get_date(project, API_KEY)

    # Now do an API search for all other repos that use showyourwork
    for page in range(1, maxpages + 1):
        req = Request(
            f"https://api.github.com/search/code?q=path:{path}+filename:{filename}&per_page=10&page={page}"
        )
        req.add_header("Accept", "application/vnd.github.v3+json")
        req.add_header("Authorization", f"token {API_KEY}")
        req.add_header("User-Agent", "request")

        # API search requests sometimes fail with `HTTP Error 403: Forbidden`
        # Not sure why! Let's try a few times before giving up.
        for j in range(maxtries):
            try:
                response = urlopen(req)
                content = response.read()
                break
            except urllib.error.HTTPError as e:
                print(f"Attempt {j+1}/{maxtries}: {e}")
                time.sleep(sleep_time)
        else:
            print("Error populating the projects page.")
            break

        # Get the dict (break if we've reached the end)
        content = json.loads(content)
        if len(content["items"]) == 0:
            break

        for res in content["items"]:
            repo = res["repository"]["full_name"]
            if (repo not in exclude_repos) and (repo not in projects):
                projects[repo] = {
                    "title": "",
                    "authors": [],
                    "field": "Uncategorized",
                    "version": get_version(repo, versions, API_KEY),
                    "commits": get_commit_count(repo, API_KEY),
                    "date": get_date(repo, API_KEY),
                }

    return projects