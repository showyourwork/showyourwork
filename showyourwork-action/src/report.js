// Imports
const core = require("@actions/core");
const shell = require("shelljs");
const { makeId, exec, getInputAsArray } = require("./utils");

// Exports
module.exports = { generateReport };

/**
 * Generate the workflow report.
 *
 */
async function generateReport() {
  core.startGroup("Generate article report");
  exec("make dag && cp src/figures/dag.pdf dag.pdf");
  core.endGroup();
  return ["dag.pdf"];
}
