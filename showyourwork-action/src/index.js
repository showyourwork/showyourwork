// Imports
const core = require("@actions/core");
const shell = require("shelljs");
const { setupConda } = require("./conda");
const { buildArticle } = require("./article");
const { publishOutput } = require("./publish");

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

  } catch (error) {

    // Exit gracefully
    core.setFailed(error.message);

  }
})();
