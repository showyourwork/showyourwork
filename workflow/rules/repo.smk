import json


localrules: repo


rule repo:
    """
    Generate repository metadata (git url, branch name, commit sha).

    """
    message:
        "Generating repo metadata..."
    output:
        posix(relpaths.temp / "repo.json"),
    priority: 99
    run:
        repo = {}
        repo["sha"] = get_repo_sha()
        repo["url"] = get_repo_url()
        repo["branch"] = get_repo_branch()
        with open(relpaths.temp / "repo.json", "w") as f:
            print(json.dumps(repo, indent=4), file=f)
