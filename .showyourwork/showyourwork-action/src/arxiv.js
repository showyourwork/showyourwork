// Imports
const core = require("@actions/core");
const shell = require("shelljs");
const { makeId, exec, getInputAsArray } = require("./utils");

// Exports
module.exports = { buildArxiv };

/**
 * Generate the workflow report.
 *
 */
async function buildArxiv() {
  if (core.getInput("upload-arxiv-artifact") == "true") {
    core.startGroup("Generate ArXiV folder");
    exec("snakemake -c1 --use-conda arxiv");
    core.endGroup();
    return ["arxiv"];
  } else {
    return [];
  }
}
