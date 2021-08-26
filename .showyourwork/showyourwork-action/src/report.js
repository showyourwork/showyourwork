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
  if (core.getInput("generate-report") == "true") {
    core.startGroup("Generate article report");
    exec("snakemake ms.pdf --report report.html");
    core.endGroup();
    return ["report.html"];
  } else {
    return [];
  }
}
