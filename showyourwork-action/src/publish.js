// Imports
const core = require("@actions/core");
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
async function publishOutput(output, report) {
  // Force-push output to `-pdf` branch
  if (core.getInput("output-branch-suffix").length > 0) {
    core.startGroup("Uploading output");
    const suffix = core.getInput("output-branch-suffix");
    const TARGET_BRANCH = `${GITHUB_BRANCH}-${suffix}`;
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
    for (const out of report) {
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
