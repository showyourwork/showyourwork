// Imports
const core = require("@actions/core");
const shell = require("shelljs");
const { formatRepo } = require("./format_repo");
const { setupConda } = require("./conda");
const { buildArticle } = require("./article");
const { generateReport } = require("./report");
const { publishOutput } = require("./publish");
const { installTeX } = require("./tex");
const { uploadTemporaries } = require("./onerror");

(async () => {
  try {
    // Exit on failure
    shell.set("-e");

    const GITHUB_SLUG = shell.env["GITHUB_REPOSITORY"];
    if (GITHUB_SLUG == "rodluger/showyourwork-template") {
      // This is a template repository -- don't do anything!
      // The workflow should be disabled by default on this
      // repo, so this is just a failsafe
    } else {
      // This is a clone of the template; let's build the paper

      // Install TeX?
      if (core.getInput("install-tex") == "true") {
        installTeX();
      }

      // Get the current showyourwork version
      const SHOWYOURWORK_VERSION = shell.exec("grep 'tag/' showyourwork/README.md").stdout.split("tag/")[1].split('"')[0];

      // Setup conda or restore from cache
      await setupConda(SHOWYOURWORK_VERSION);

      // Build the article
      var output = await buildArticle(SHOWYOURWORK_VERSION);

      // Generate the report
      var report = await generateReport();

      // Publish the article output
      await publishOutput(output, report);

      // Format repository if it's a fresh fork
      formatRepo();
    }
  } catch (error) {
    // Exit gracefully
    await uploadTemporaries();
    core.setFailed(error.message);
  }
})();
