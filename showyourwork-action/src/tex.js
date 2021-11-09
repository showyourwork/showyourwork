// Imports
const { makeId, exec, getInputAsArray } = require("./utils");

// Exports
module.exports = { installTeX };

/**
 * Setup a minimal TeX distribution.
 *
 */
function installTeX() {

  // Download and setup TinyTex
  exec(`wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh" | sh && ~/bin/latex --version`, "Download and install TinyTex");

  // Get custom package names
  const packages_array = getInputAsArray("tex-packages");
  var packages = '';
  for (var i in packages_array) {
      packages += ' ' + packages_array[i];
  }

  // Install custom packages
  try {
    exec(`sudo ~/bin/tlmgr install ${packages}`);
  } catch (error) {
    // Update tlmgr
    exec(`sudo ~/bin/tlmgr update --self --all`);
    exec(`sudo ~/bin/tlmgr path add`);
    exec(`sudo ~/bin/fmtutil-sys -all`);
    exec(`sudo ~/bin/tlmgr install ${packages}`);
  }
}
