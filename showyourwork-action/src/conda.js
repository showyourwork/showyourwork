// Imports
const core = require("@actions/core");
const cache = require("@actions/cache");
const shell = require("shelljs");
const constants = require("./constants.js");
const { exec } = require("./utils");

// Exports
module.exports = { setupConda };

// Cache settings
const CONDA_CACHE_NUMBER = core.getInput("conda-cache-number");
const RUNNER_OS = shell.env["RUNNER_OS"];
const conda_key = `conda-${constants.conda_cache_version}-${RUNNER_OS}-${CONDA_CACHE_NUMBER}`;
const conda_restoreKeys = [];
const conda_paths = ["~/.conda", "~/.condarc", "~/conda_pkgs_dir", "sywenvs"];

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
  if (!shell.test("-d", "./sywenvs")) {
    exec(
      ". ~/.conda/etc/profile.d/conda.sh && conda create -y -p ./sywenvs",
      "Create environment"
    );
    exec(
      "make install_deps",
      "Install dependencies"
    );
  }

  // Save conda cache (failure OK)
  try {
    core.startGroup("Update conda cache");
    const conda_cacheId = await cache.saveCache(conda_paths, conda_key);
    core.endGroup();
  } catch (error) {
    shell.echo(`WARNING: ${error.message}`);
  }
}
