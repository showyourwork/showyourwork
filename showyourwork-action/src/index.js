// Imports
const core = require("@actions/core");
const shell = require("shelljs");
const { setupConda } = require("./conda");
const { buildArticle } = require("./article");
const { buildTarball } = require("./arxiv");
const { publishOutput } = require("./publish");
const { publishLogs } = require("./logs");

(async () => {
  try {
    
    // Exit on failure
    shell.set("-e");

    // Setup conda or restore from cache
    await setupConda();

    // Build the article
    await buildArticle();

    // Build arxiv tarball
    if (core.getInput("arxiv-tarball") == "true") {
      await buildTarball();
    }

    // Publish the article output
    await publishOutput();

    // Publish the logs
    await publishLogs();

  } catch (error) {

    // Publish the logs
    try {
      await publishLogs();
    } catch (error) {
      core.error("Unable to upload the build logs.");
      core.setFailed(error.message);
    }

    // Exit gracefully
    core.setFailed(error.message);

  }
})();
