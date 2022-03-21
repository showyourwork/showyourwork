// Imports
const core = require("@actions/core");
const shell = require("shelljs");

// Exports
module.exports = { makeId, exec, getInputAsArray };

/**
 * Generate a random hash.
 *
 * See https://stackoverflow.com/a/1349426
 *
 * and
 *
 * https://github.com/actions/cache/issues/432#issuecomment-740376179
 */
function makeId(length) {
  var result = "";
  var characters = "abcdefghijklmnopqrstuvwxyz0123456789";
  var charactersLength = characters.length;
  for (var i = 0; i < length; i++) {
    result += characters.charAt(Math.floor(Math.random() * charactersLength));
  }
  return result;
}

/**
 * Simple wrapper to execute a bash command, optionally in a log file group.
 *
 */
function exec_wrapper(cmd, group) {
  if (typeof group !== "undefined") {
    core.startGroup(group);
  }
  const result = shell.exec(cmd, { shell: "/bin/bash" });
  if (result.code != 0) {
    shell.echo(`Error: ${cmd}`);
    shell.exit(1);
  }
  if (typeof group !== "undefined") {
    core.endGroup();
  }
  return result;
}

/**
 * Execute a bash command in conda environment `./sywenvs`, if present.
 *
 */
function exec(cmd, group) {
  if (
    shell.test("-f", "~/.conda/etc/profile.d/conda.sh") &&
    shell.test("-d", "./sywenvs")
  ) {
    return exec_wrapper(
      `. ~/.conda/etc/profile.d/conda.sh && conda activate ./sywenvs && ${cmd}`,
      group
    );
  } else {
    return exec_wrapper(cmd, group);
  }
}

/**
 * Get a YAML input as an array.
 *
 */
function getInputAsArray(name) {
  return core
    .getInput(name)
    .split("\n")
    .map((s) => s.trim())
    .filter((x) => x !== "");
}
