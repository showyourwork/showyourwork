// Imports
const core = require("@actions/core");
const artifact = require("@actions/artifact");
const shell = require("shelljs");

// Exports
module.exports = { publishOutput };

// Get repo info
const GITHUB_SLUG = shell.env["GITHUB_REPOSITORY"];
const GITHUB_BRANCH = shell
  .exec("git rev-parse --abbrev-ref HEAD")
  .replace(/(\r\n|\n|\r)/gm, "");
const GITHUB_TOKEN = core.getInput("github-token");
const GITHUB_WORKSPACE = shell.env["GITHUB_WORKSPACE"];

/**
 * Publish the article output.
 *
 */
async function publishOutput(output, arxiv, report) {
  // Upload artifact
  if (core.getInput("upload-arxiv-artifact") == "true") {
    core.startGroup("Upload arxiv artifact");
    const artifactClient = artifact.create();
    const artifactName = "arxiv";
    const rootDirectory = ".";
    const options = {
      continueOnError: false,
    };
    const uploadResponse = await artifactClient.uploadArtifact(
      artifactName,
      arxiv,
      rootDirectory,
      options
    );
    core.endGroup();
  }

  // Force-push to `-pdf` branch
  if (core.getInput("force-push") == "true") {
    core.startGroup("Uploading output");
    const TARGET_BRANCH = `${GITHUB_BRANCH}-pdf`;
    const TARGET_DIRECTORY = shell
      .exec("mktemp -d")
      .replace(/(\r\n|\n|\r)/gm, "");
    shell.cp("-R", ".", `${TARGET_DIRECTORY}`);
    shell.cd(`${TARGET_DIRECTORY}`);
    shell.exec(`git checkout --orphan ${TARGET_BRANCH}`);
    var silentState = shell.config.silent;
    shell.config.silent = true;
    shell.exec("git rm --cached -rf .");
    shell.config.silent = silentState;
    for (const out of output) {
      shell.exec(`git add -f ${out}`);
    }
    shell.exec(
      "git -c user.name='showyourwork' -c user.email='showyourwork' " +
        "commit -m 'force-push article output'"
    );
    shell.exec(
      "git push --force " +
        `https://x-access-token:${GITHUB_TOKEN}@github.com/${GITHUB_SLUG} ` +
        `${TARGET_BRANCH}`
    );
    shell.cd(GITHUB_WORKSPACE);
    core.endGroup();
  }

  // Upload the report to GitHub Pages
  if (report.length > 0) {
    core.startGroup("Uploading report");
    const TARGET_BRANCH = core.getInput("gh-pages-branch");
    if (TARGET_BRANCH == GITHUB_BRANCH) {
      core.setFailed("GitHub Pages branch can't be the current branch!");
    }
    const TARGET_DIRECTORY = shell
      .exec("mktemp -d")
      .replace(/(\r\n|\n|\r)/gm, "");
    shell.cd(`${TARGET_DIRECTORY}`);
    // Attempt to pull from the remote, if the GitHub Pages branch exists
    shell.set("+e");
    const result = shell.exec(
      `git clone https://github.com/${GITHUB_SLUG} --branch ${TARGET_BRANCH} --single-branch .`
    );
    shell.set("-e");
    if (result.code == 0) {
      // Branch exists
      for (const out of report) {
        shell.exec(`cp -r ${GITHUB_WORKSPACE}/${out} .`);
        shell.exec(`git add -f ${out}`);
      }
      shell.exec(
        "git -c user.name='showyourwork' -c user.email='showyourwork' " +
          "commit -m 'update article report'"
      );
      shell.exec(
        "git push " +
          `https://x-access-token:${GITHUB_TOKEN}@github.com/${GITHUB_SLUG} ` +
          `${TARGET_BRANCH}`
      );
    } else {
      // Branch doesn't exist (create it & force push)
      shell.exec("git init");
      shell.exec(`git checkout --orphan ${TARGET_BRANCH}`);
      for (const out of report) {
        shell.exec(`cp -r ${GITHUB_WORKSPACE}/${out} .`);
        shell.exec(`git add -f ${out}`);
      }
      shell.exec(
        "git -c user.name='showyourwork' -c user.email='showyourwork' " +
          "commit -m 'create gh-pages branch and upload article report'"
      );
      shell.exec(
        "git push --force " +
          `https://x-access-token:${GITHUB_TOKEN}@github.com/${GITHUB_SLUG} ` +
          `${TARGET_BRANCH}`
      );
    }
    shell.cd(GITHUB_WORKSPACE);
    core.endGroup();
  }
}
