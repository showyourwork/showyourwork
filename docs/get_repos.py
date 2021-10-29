import urllib
from urllib.request import Request, urlopen
import json
import os
import sys


def get_repos(
    filename="showyourwork.yml",
    path=".github/workflows",
    maxpages=10,
    exclude_repos=[
        "rodluger/showyourwork-template",
        "rodluger/showyourwork-sandbox",
        "rodluger/showyourwork-example-dev",
    ],
):
    """
    Return a list of repos that use showyourwork.

    """
    # Get curated list of projects that use showyourwork
    with open("projects.json", "r") as f:
        curated = json.load(f)

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
    projects = {}
    for page in range(1, maxpages + 1):
        req = Request(
            f"https://api.github.com/search/code?q=path:{path}+filename:{filename}&per_page=100&page={page}"
        )
        req.add_header("Accept", "application/vnd.github.v3+json")
        req.add_header("Authorization", f"token {API_KEY}")
        content = urlopen(req).read()
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

                # Get the topic and description
                for topic in curated:
                    if repo in curated[topic]:
                        description = curated[topic][repo]
                        break
                else:
                    topic = "Uncategorized"
                    description = ""

                # Assemble the metadata
                if topic not in projects:
                    projects[topic] = []
                projects[topic].append(
                    {
                        "name": repo,
                        "description": description,
                        "version": version,
                        "date": date,
                    }
                )

    # Sort by the timestamps
    for topic in projects:
        projects[topic] = sorted(projects[topic], key=lambda item: item["date"])[::-1]

    breakpoint()

    return projects
