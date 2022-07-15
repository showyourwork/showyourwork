from helpers import TemporaryShowyourworkRepository
from showyourwork import gitapi
from showyourwork.subproc import get_stdout, parse_request
import pytest
import requests
import json
import asyncio
import re

@pytest.mark.xfail("Error in latexdiff: You can't use `/raise' in internal vertical mode")
class TestPullRequest(TemporaryShowyourworkRepository):
    """
    Test pull request behavior, including the generation of the diff and
    the automated PR message w/ the link to the PDF.

    """

    def build_local(self, pre=""):
        """No need to build the article locally."""
        pass

    def git_commit(self):
        """Nothing to commit initially."""
        pass

    @pytest.mark.asyncio_cooperative
    async def run_github_action(self):
        """Make changes in a new branch, issue a PR, and inspect the rendered diff."""
        # Push the first commit to the main branch (w/out running the action)
        get_stdout(
            "git commit --amend -m '[skip ci] auto commit from showyourwork tests'",
            shell=True,
            cwd=self.cwd,
        )
        get_stdout(
            "git push --force https://x-access-token:"
            f"{gitapi.get_access_token()}"
            f"@github.com/showyourwork/{self.repo} main",
            shell=True,
            cwd=self.cwd,
            secrets=[gitapi.get_access_token()],
        )

        # Create a new branch
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

        # Wait for the action to complete & check its status
        await asyncio.sleep(self.action_wait)
        status = "unknown"
        for n in range(self.action_max_tries):

            # Note that we're querying `workflow_run` actions, specifically
            # to catch the `process-pull-request` action and check that the bot
            # has posted a comment to the PR with the link to the PDF and the diff.
            (status, conclusion, url,) = gitapi.get_latest_workflow_run_status(
                self.repo, org="showyourwork", event="workflow_run"
            )
            if status == "completed":
                if conclusion == "success":
                    print(f"[{self.repo}] Workflow completed successfully.")
                    return
                else:
                    raise Exception(
                        "[{self.repo}] GitHub Actions workflow terminated "
                        f"with status {conclusion}.\n"
                        f"For details, see {url}"
                    )
            elif n < self.action_max_tries - 1:
                print(
                    f"[{self.repo}] Waiting {self.action_interval} seconds for "
                    f"workflow to finish ({n+2}/{self.action_max_tries})..."
                )
                await asyncio.sleep(self.action_interval)
        else:
            raise Exception(
                "[{self.repo}] GitHub Actions workflow timed out.\n"
                f"For details, see {url}"
            )

        # Get the bot comment
        url = f"https://api.github.com/repos/showyourwork/{self.repo}/issues/1/comments"
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
                        "\[PDF with highlighted changes\]\((.*?)\)", comment
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
        