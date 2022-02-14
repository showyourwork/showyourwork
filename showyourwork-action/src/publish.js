// Imports
const core = require("@actions/core");
const artifact = require("@actions/artifact");
const shell = require("shelljs");

// Exports
module.exports = { publishOutput };

// Get repo info
const GITHUB_SLUG = shell.env["GITHUB_REPOSITORY"];
const GITHUB_BRANCH = shell
  .exec("echo ${GITHUB_REF##*/}", {silent: true})
  .replace(/(\r\n|\n|\r)/gm, "");
const GITHUB_TOKEN = core.getInput("github-token");
const GITHUB_WORKSPACE = shell.env["GITHUB_WORKSPACE"];
const PULL_REQUEST = shell.env["GITHUB_BASE_REF"].length > 0;

/**
 * Publish the article output.
 *
 */
async function publishOutput() {

  // Infer the manuscript name
  const config = require(`${GITHUB_WORKSPACE}/.showyourwork/config.json`);
  const output = [config["ms_pdf"]];

  // Upload the arxiv tarball?
  if (core.getInput("arxiv-tarball") == "true") {
    output.push("arxiv.tar.gz");
  }

  if (PULL_REQUEST) {
    // Upload an artifact
    const artifactClient = artifact.create();
    const uploadResponse = await artifactClient.uploadArtifact(
      "showyourwork-output", 
      output, 
      ".", 
      {
        continueOnError: false
      }
    );
  } else {
    // Force-push output to a separate branch
    if ((core.getInput("output-branch-suffix").length > 0) && (GITHUB_BRANCH.length > 0)) {
      core.startGroup("Uploading output");
      const suffix = core.getInput("output-branch-suffix");
      const TARGET_BRANCH = `${GITHUB_BRANCH}-${suffix}`;
      const TARGET_DIRECTORY = shell
        .exec("mktemp -d")
        .replace(/(\r\n|\n|\r)/gm, "");
      shell.cp("-R", ".", `${TARGET_DIRECTORY}`);
      shell.cd(`${TARGET_DIRECTORY}`);
      shell.exec(`git checkout --orphan ${TARGET_BRANCH}`);
      shell.exec("git rm --cached -rf .", {silent: true});
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
  }
}
