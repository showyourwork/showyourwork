// Imports
const core = require("@actions/core");
const cache = require("@actions/cache");
const shell = require("shelljs");
const { makeId, exec, getInputAsArray } = require("./utils");

// Exports
module.exports = { setupConda };

// Cache settings
const CONDA_CACHE_NUMBER = core.getInput("conda-cache-number");
const RUNNER_OS = shell.env["RUNNER_OS"];
const conda_key = `debug3-conda-${RUNNER_OS}-${CONDA_CACHE_NUMBER}`; // DEBUG
const conda_restoreKeys = [];
const conda_paths = ["~/.conda", "~/.condarc", "~/conda_pkgs_dir", "envs"];

/**
 * Setup a conda distribution or restore it from cache.
 *
 */
async function setupConda() {
  // Restore conda cache
  core.startGroup("Restore conda cache");
  const conda_cacheKey = await cache.restoreCache(
    conda_paths,
    conda_key,
    conda_restoreKeys
  );
  core.endGroup();

  // Download and setup conda
  if (!shell.test("-d", "~/.conda")) {
    exec(
      "wget --no-verbose https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ./conda.sh", 
      "Download conda"
    );
    exec("bash ./conda.sh -b -p ~/.conda && rm -f ./conda.sh", "Install conda");
    exec(
      ". ~/.conda/etc/profile.d/conda.sh && " +
      "conda config --add pkgs_dirs ~/conda_pkgs_dir"
    );
  }

  // Display some info
  exec(". ~/.conda/etc/profile.d/conda.sh && conda info", "Conda info");

  // Create environment & install snakemake
  if (!shell.test("-d", "./envs")) {
    exec(
      ". ~/.conda/etc/profile.d/conda.sh && conda create -y -p ./envs",
      "Create environment"
    );
    exec(
      "make snakemake_setup",
      "Install snakemake"
    );
  }

  // Save conda cache (failure OK)
  try {
    core.startGroup("Update conda cache");
    const conda_cacheId = await cache.saveCache(conda_paths, conda_key);
    core.endGroup();
  } catch (error) {
    core.warning(error.message);
  }
}
