// Imports
const core = require("@actions/core");
const artifact = require("@actions/artifact");
const shell = require("shelljs");
const glob = require("glob");

// Exports
module.exports = { uploadTemporaries };

/**
 * Upload build logs and temporaries for debugging.
 *
 */
async function uploadTemporaries() {

    // Maximum size of individual files; files bigger
    // than this will be skipped.
    const maxSizeInKB = 5000;

    // Assemble a list of folders/files to upload
    const patterns = [
      "*.*",
      "src/*.*"
    ];
    var files = [];
    patterns.forEach(function (pattern) {
      glob(pattern, (err, fs) => {
        fs.forEach(function (f) {
          try {
            var sz = parseInt(shell.exec(`du -k ${f} | cut -f1`, {silent:true}).stdout.replace('\n', ''));
            if (sz < maxSizeInKB) {
              files.push(f);
            }
            
            // DEBUG
            core.info(f);
            core.info(sz);
            core.info("<FILES>");
            core.info(files);
            core.info("</FILES>");

          } catch (error) {
            core.warning(error.message);
          }
        })
      })
    });

    // DEBUG
    core.info("<FILES>");
    core.info(files);
    core.info("</FILES>");

    // DEBUG
    shelljs.exec("more showyourwork/Makefile");

    // Upload the artifact
    const artifactClient = artifact.create();
    const uploadResponse = await artifactClient.uploadArtifact(
      "showyourwork-output", 
      files, 
      ".", 
      {
        continueOnError: true
      }
    );
   
}
