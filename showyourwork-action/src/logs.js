// Imports
const core = require("@actions/core");
const artifact = require("@actions/artifact");
const shell = require("shelljs");

// Exports
module.exports = { publishLogs };


/**
 * Publish an artifact containing all the showyourwork logs and temp files.
 *
 */
async function publishLogs() {

    core.startGroup("Uploading logs");

    // Collect all files to upload
    var files = shell.exec(
        "find .showyourwork " + 
        "-not -type d " + 
        "-not -path '.showyourwork/cache/**' " + 
        "-not -path '.showyourwork/zenodo/**' " + 
        "-not -path '.showyourwork/zenodo_sandbox/**' " + 
        "-print", 
        {silent: true}
    );
    files = files.split("\n").filter(n => n);

    // Upload the artifact
    const artifactClient = artifact.create();
    const uploadResponse = await artifactClient.uploadArtifact(
        "showyourwork-logs", 
        files, 
        ".showyourwork", 
        {
            continueOnError: false
        }
    );
    core.endGroup();

    core.startGroup("Repository tree");
    shell.exec("tree -DaL 1");
    shell.exec("tree -Da src");
    core.endGroup();
}
