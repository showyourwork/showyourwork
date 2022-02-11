// Imports
const core = require("@actions/core");
const shell = require("shelljs");
const { setupConda } = require("./conda");
const { buildArticle } = require("./article");
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

    // Publish the article output
    await publishOutput();

    // Publish the logs
    await publishLogs();

  } catch (error) {

    // Publish the logs
    await publishLogs();

    // Exit gracefully
    core.setFailed(error.message);

  }
})();
