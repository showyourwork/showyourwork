// Imports
const core = require("@actions/core");
const cache = require("@actions/cache");
const {makeId, exec} = require("./utils");


// Exports
module.exports = {buildArticle};


/**
 * Build the article.
 *
 */
 function buildArticle() {

    // TODO: Restore article cache

    // Outputs
    var output = [];

    // Build the article
    core.startGroup("Build article");
    if (core.getInput("verbose") == "true") {
        exec("snakemake -c1 --use-conda ms.pdf --verbose");
    } else {
        exec("snakemake -c1 --use-conda ms.pdf");
    }
    output.push("ms.pdf");
    core.endGroup();

    // Generate DAG?
    if (core.getInput("generate-dag") == "true") {
      core.startGroup("Generate article DAG");
      exec("snakemake ms.pdf --dag | dot -Tpdf > dag.pdf");
      output.push("dag.pdf");
      core.endGroup();
    }

    // TODO: Save article cache

    return output;

 }
