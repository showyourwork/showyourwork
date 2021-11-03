// Imports
const core = require("@actions/core");
const shell = require("shelljs");

// Exports
module.exports = { installTeX };

/**
 * Setup a minimal TeX distribution.
 *
 */
function installTeX() {

  // Download and setup TinyTex
  exec(`wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh" | sh && ~/bin/latex --version`, "Download and install TinyTex");

}
