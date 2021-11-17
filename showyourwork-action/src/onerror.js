// Imports
const core = require("@actions/core");
const artifact = require("@actions/artifact");
const shell = require("shelljs");

// Exports
module.exports = { uploadTemporaries };

/**
 * Upload build logs and temporaries for debugging.
 *
 */
async function uploadTemporaries() {

    // TODO: Assemble a list of files to upload

    // Upload an artifact
    const files = {};
    const artifactClient = artifact.create();
    const uploadResponse = await artifactClient.uploadArtifact(
      "showyourwork-output", 
      files, 
      ".", 
      {
        continueOnError: false
      }
    );
}
