// Imports
const artifact = require("@actions/artifact");
const shell = require("shelljs");

// Exports
module.exports = { publishLogs };


/**
 * Publish an artifact containing all the showyourwork logs and temp files.
 *
 */
async function publishLogs() {
    // Upload an artifact
    shell.exec("tar -czvf logs.tar.gz .showyourwork");
    const artifactClient = artifact.create();
    const uploadResponse = await artifactClient.uploadArtifact(
        "showyourwork-logs", 
        ["logs.tar.gz"], 
        ".", 
        {
            continueOnError: false
        }
    );
}
