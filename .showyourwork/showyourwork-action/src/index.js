// Imports
const core = require("@actions/core");
const shell = require("shelljs");
const {formatRepo} = require("./format_repo");

(async () => {

  try {

    // Exit on failure
    shell.set("-e");

    // Format repository if it's a fresh fork
    formatRepo();

  } catch (error) {

    // Exit gracefully
    core.setFailed(error.message);

  }

})();
