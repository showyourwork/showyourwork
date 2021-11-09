// Imports
const core = require("@actions/core");
const cache = require("@actions/cache");
const shell = require("shelljs");
const { makeId, exec, getInputAsArray } = require("./utils");

// Exports
module.exports = { buildArticle };

/**
 * Build the article.
 *
 */
async function buildArticle(ARTICLE_CACHE_NUMBER = null) {
  // Article cache settings. We're caching pretty much everything
  // in the repo, but overriding it with any files that changed since
  // the commit at which we cached everything
  // Note that the GITHUB_REF (branch) is part of the cache key
  // so we don't mix up the caches for different branches!
  const ACTION_PATH = core.getInput("action-path");
  if (ARTICLE_CACHE_NUMBER == null)
    ARTICLE_CACHE_NUMBER = core.getInput("article-cache-number");
  const RUNNER_OS = shell.env["RUNNER_OS"];
  const GITHUB_REF = shell.env["GITHUB_REF"];
  const randomId = makeId(8);
  const article_key = `article-${RUNNER_OS}-${GITHUB_REF}-${ARTICLE_CACHE_NUMBER}-${randomId}`;
  const article_restoreKeys = [
    `article-${RUNNER_OS}-${GITHUB_REF}-${ARTICLE_CACHE_NUMBER}`,
  ];
  const article_paths = [
    ".snakemake",
    ".showyourwork",
    ".last-commit",
    "environment.yml",
    "ms.pdf",
    "src",
  ];

  // Restore the article cache
  core.startGroup("Restore article cache");
  const article_cacheKey = await cache.restoreCache(
    article_paths,
    article_key,
    article_restoreKeys
  );
  exec(`rm -f .showyourwork/repo.json`); // Always re-generate this!
  exec(`python ${ACTION_PATH}/src/cache.py --restore`);
  core.endGroup();

  // Outputs
  var output = [];

  // Build the article
  core.startGroup("Build article");
  if (core.getInput("verbose") == "true") {
    exec("make ms.pdf OPTIONS='-c1 --verbose --reason --notemp'");
  } else {
    exec("make ms.pdf OPTIONS='-c1 --reason --notemp'");
  }
  output.push("ms.pdf");
  core.endGroup();

  // Build arxiv tarball
  if (core.getInput("arxiv-tarball") == "true") {
    core.startGroup("Build ArXiV tarball");
    if (core.getInput("verbose") == "true") {
      exec(
        `make arxiv.tar.gz OPTIONS='-c1 --verbose --reason --notemp'`
      );
    } else {
      exec(
        `make arxiv.tar.gz OPTIONS='-c1 --reason --notemp'`
      );
    }
    output.push("arxiv.tar.gz");
    core.endGroup();
  }

  

  // Save article cache
  core.startGroup("Update article cache");
  exec("make remove_zenodo_datasets");
  exec(`python ${ACTION_PATH}/src/cache.py --update`);
  const article_cacheId = await cache.saveCache(article_paths, article_key);
  core.endGroup();

  return output;
}
