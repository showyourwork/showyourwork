// Imports
const core = require("@actions/core");
const shell = require("shelljs");
const {formatRepo} = require("./format_repo");
const {setupConda} = require("./conda");
const {buildArticle} = require("./article");
const {publishOutput} = require("./publish");


(async () => {

  try {

    // Exit on failure
    shell.set("-e");

    // Format repository if it's a fresh fork
    formatRepo();

    // Setup conda or restore from cache
    await setupConda();

    // Build the article
    output = buildArticle();

    // Publish the article output
    publishOutput(output);

  } catch (error) {

    // Exit gracefully
    core.setFailed(error.message);

  }

})();
