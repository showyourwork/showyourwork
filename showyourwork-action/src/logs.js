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

    // Save a snapshot of the repository tree
    shell.exec("tree -Dah --timefmt='%d/%m/%Y %H:%M:%S' -I 'showyourwork|.snakemake|.git|__pycache__|sywenvs' > .showyourwork/tree.txt");

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

}
