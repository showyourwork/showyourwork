import asyncio
import json
import re

import pytest
import requests
from helpers import TemporaryShowyourworkRepository

from showyourwork import gitapi
from showyourwork.subproc import get_stdout, parse_request

pytestmark = pytest.mark.remote


class TestPullRequests(TemporaryShowyourworkRepository):
    """
    Test pull request behavior, including the generation of the diff using
    ``latex-diff`` and the automated PR message with the link to the PDF.

    """

    @pytest.mark.asyncio_cooperative
    async def run_github_action(self):
        """Make changes in a new branch, issue a PR, and inspect the rendered diff."""
        # Push the initial commit to main
        print(f"[{self.repo}] Pushing to `showyourwork/{self.repo}@main`...")
        get_stdout(
            "git push --force https://x-access-token:"
            f"{gitapi.get_access_token()}"
            f"@github.com/showyourwork/{self.repo} main",
            shell=True,
            cwd=self.cwd,
            secrets=[gitapi.get_access_token()],
        )
        # Avoid pushing two commits simultaneously, since this was causing
        # the Actions to occasionally not run on the PR.
        await asyncio.sleep(30)
        head_sha = get_stdout(
            "git rev-parse HEAD", shell=True, cwd=self.cwd
        ).replace("\n", "")

        # Create a new branch
        print(
            f"[{self.repo}] Creating branch `small-change` with some changes..."
        )
        get_stdout("git checkout -b small-change", shell=True, cwd=self.cwd)

        # Make a few small changes, deletions, and insertions.
        with open(self.cwd / "src" / "tex" / "ms.tex", "r") as f:
            contents = f.read()
        contents = contents.replace(
            "An open source scientific article",
            "An awesome open source scientific article",
        )
        contents = contents.replace(
            "Ut purus elit, vestibulum ut, placerat",
            "Placerat",
        )
        contents = contents.replace("Lorem ipsum", "Merol muspi")
        with open(self.cwd / "src" / "tex" / "ms.tex", "w") as f:
            f.write(contents)

        # Add, commit, and push to the new branch
        print(
            f"[{self.repo}] Pushing to `showyourwork/{self.repo}@small-change`..."
        )
        get_stdout("git add .", shell=True, cwd=self.cwd)
        get_stdout(
            "git -c user.name='gh-actions' -c user.email='gh-actions' "
            "commit -q -m 'auto commit from showyourwork tests'",
            shell=True,
            cwd=self.cwd,
        )
        get_stdout(
            "git push --force https://x-access-token:"
            f"{gitapi.get_access_token()}"
            f"@github.com/showyourwork/{self.repo} small-change",
            shell=True,
            cwd=self.cwd,
            secrets=[gitapi.get_access_token()],
        )

        # Create the PR for the main branch
        print(f"[{self.repo}] Creating PR for small-change -> main...")
        url = f"https://api.github.com/repos/showyourwork/{self.repo}/pulls"
        data = parse_request(
            requests.post(
                url,
                headers={
                    "Accept": "application/vnd.github.v3+json",
                    "Authorization": f"token {gitapi.get_access_token()}",
                },
                data=json.dumps(
                    {
                        "title": "Pull request test",
                        "body": "A test of the diff feature of showyourwork pull requests.",
                        "head": "small-change",
                        "base": "main",
                    }
                ),
            )
        )
        pr_number = data["number"]

        # Wait for the action to complete & check its status
        print(
            f"[{self.repo}] Waiting {self.action_wait} seconds for PR workflow "
            f"to finish (1/{self.action_max_tries})..."
        )
        await asyncio.sleep(self.action_wait)
        status = "unknown"
        for n in range(self.action_max_tries):

            # Note that we're querying `workflow_run` actions, specifically
            # to catch the `process-pull-request` action and check that the bot
            # has posted a comment to the PR with the link to the PDF and the diff.
            (status, conclusion, url,) = gitapi.get_workflow_run_status(
                self.repo,
                org="showyourwork",
                q={
                    "event": "workflow_run",
                    "head_sha": head_sha,
                },
            )
            if status == "completed":
                if conclusion == "success":
                    print(f"[{self.repo}] Workflow completed successfully.")
                    break
                else:
                    raise Exception(
                        "[{self.repo}] GitHub Actions workflow terminated "
                        f"with status {conclusion}.\n"
                        f"For details, see {url}"
                    )
            elif n < self.action_max_tries - 1:
                print(
                    f"[{self.repo}] Waiting {self.action_interval} seconds for "
                    f"PR workflow to finish ({n+2}/{self.action_max_tries})..."
                )
                await asyncio.sleep(self.action_interval)
        else:
            raise Exception(
                "[{self.repo}] GitHub Actions workflow timed out.\n"
                f"For details, see {url}"
            )

        # Get the bot comment
        url = f"https://api.github.com/repos/showyourwork/{self.repo}/issues/{pr_number}/comments"
        data = parse_request(
            requests.get(
                url,
                headers={
                    "Accept": "application/vnd.github.v3+json",
                    "Authorization": f"token {gitapi.get_access_token()}",
                },
            )
        )
        for comment in data:
            if comment["user"]["login"] == "github-actions[bot]":
                comment = comment["body"]
                try:
                    diff_url = re.search(
                        r"\[PDF with highlighted changes\]\((.*?)\)", comment
                    ).groups()[0]
                except:
                    raise Exception(
                        "Bot did not post the link to the PDF diff."
                    )
                else:
                    break
        else:
            raise Exception("Cannot find bot comment on PR.")

        # Download the PDF diff
        response = requests.get(diff_url)
        if response.status_code > 204:
            raise Exception("Diff PDF was not pushed to the -pdf branch.")
        with open(self.cwd / "diff.pdf", "wb") as f:
            f.write(response.content)

        # TODO: Inspect the diff
        pass

        # Close the PR
        print(f"[{self.repo}] Closing the PR...")
        url = f"https://api.github.com/repos/showyourwork/{self.repo}/pulls/{pr_number}"
        data = parse_request(
            requests.patch(
                url,
                headers={
                    "Accept": "application/vnd.github.v3+json",
                    "Authorization": f"token {gitapi.get_access_token()}",
                },
                data=json.dumps({"state": "closed"}),
            )
        )
