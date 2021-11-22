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
  exec(
    "sudo apt-get install graphviz",
    "Install graphviz"
  );
  core.startGroup("Generate article report");
  exec("snakemake ms.pdf --dag | dot -Tpdf > dag.pdf");
  core.endGroup();
  return ["dag.pdf"];
}
