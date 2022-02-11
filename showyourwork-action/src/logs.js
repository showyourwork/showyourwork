// Imports
const artifact = require("@actions/artifact");

// Exports
module.exports = { publishLogs };


/**
 * Publish an artifact containing all the showyourwork logs and temp files.
 *
 */
async function publishLogs() {
    // Upload an artifact
    const artifactClient = artifact.create();
    const uploadResponse = await artifactClient.uploadArtifact(
        "showyourwork-output", 
        [".showyourwork"], 
        ".", 
        {
            continueOnError: false
        }
    );
}
