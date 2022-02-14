// Imports
const core = require("@actions/core");
const { exec } = require("./utils");

// Exports
module.exports = { buildTarball };

/**
 * Build the arXiv tarball.
 *
 */
async function buildTarball() {

  core.startGroup("Build arxiv tarball");
  exec("make arxiv");
  core.endGroup();

}
