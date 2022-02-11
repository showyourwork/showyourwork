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
    const files = shell.exec("find .showyourwork -not -type d -print", {silent: true}).split("\n").filter(n => n);
    const artifactClient = artifact.create();
    const uploadResponse = await artifactClient.uploadArtifact(
        "showyourwork-logs", 
        files, 
        ".", 
        {
            continueOnError: false
        }
    );
}
