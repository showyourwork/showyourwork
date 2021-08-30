// Imports
const core = require("@actions/core");
const shell = require("shelljs");

// Exports
module.exports = { runTests };

// Get repo info
const GITHUB_SLUG = shell.env["GITHUB_REPOSITORY"];
const GITHUB_BRANCH = shell
  .exec("git rev-parse --abbrev-ref HEAD")
  .replace(/(\r\n|\n|\r)/gm, "");
const GITHUB_TOKEN = core.getInput("github-token");
const GITHUB_WORKSPACE = shell.env["GITHUB_WORKSPACE"];

const ACCESS_TOKEN = core.getInput("access-token");
core.setSecret(ACCESS_TOKEN);

/**
 * Run tests.
 *
 */
async function runTests(output, report) {
  const TEST_SLUG = `${GITHUB_SLUG}-test`;
  core.startGroup(`Creating ${TEST_SLUG}`);
  const TARGET_DIRECTORY = shell
    .exec("mktemp -d")
    .replace(/(\r\n|\n|\r)/gm, "");
  shell.cp("-R", ".", `${TARGET_DIRECTORY}`);
  shell.cd(`${TARGET_DIRECTORY}`);
  shell.exec("rm -rf .git");
  shell.exec("git init");
  shell.exec(`git checkout --orphan ${TARGET_BRANCH}`);
  var silentState = shell.config.silent;
  shell.config.silent = true;
  shell.exec("git add .");
  shell.exec(
    "git -c user.name='gh-actions' -c user.email='gh-actions' commit -m 'test commit'"
  );
  shell.config.silent = silentState;
  shell.exec(
    "git push --force " +
      `https://x-access-token:${ACCESS_TOKEN}@github.com/${TEST_SLUG} main`
  );
  shell.cd(GITHUB_WORKSPACE);
  core.endGroup();
}
