// Imports
const core = require("@actions/core");
const cache = require("@actions/cache");
const shell = require("shelljs");
const { makeId, exec, getInputAsArray } = require("./utils");

// Exports
module.exports = { buildArticle };

// Article cache settings. We're caching pretty much everything
// in the repo, but overriding it with any files that changed since
// the commit at which we cached everything
const ACTION_PATH = core.getInput("action-path");
const ARTICLE_CACHE_NUMBER = core.getInput("article-cache-number");
const RUNNER_OS = shell.env["RUNNER_OS"];
const randomId = makeId(8);
const article_key = `article-${RUNNER_OS}-${ARTICLE_CACHE_NUMBER}-${randomId}`;
const article_restoreKeys = [`article-${RUNNER_OS}-${ARTICLE_CACHE_NUMBER}`];
const article_paths = getInputAsArray("article-cache-paths").concat([
  ".snakemake",
  ".showyourwork/tmp",
  ".showyourwork/resources",
  "environment.yml",
  "ms.pdf",
  "figures",
  "data",
  "tex",
]);

/**
 * Build the article.
 *
 */
async function buildArticle() {
  // Restore the article cache
  core.startGroup("Restore article cache");
  const article_cacheKey = await cache.restoreCache(
    article_paths,
    article_key,
    article_restoreKeys
  );
  exec(`python ${ACTION_PATH}/src/cache.py --restore`);
  core.endGroup();

  // Outputs
  var output = [];

  // Build the article
  core.startGroup("Build article");
  if (core.getInput("verbose") == "true") {
    exec("snakemake -c1 --use-conda ms.pdf --verbose --reason");
  } else {
    exec("snakemake -c1 --use-conda ms.pdf --reason");
  }
  output.push("ms.pdf");
  core.endGroup();

  // Save article cache (failure OK)
  try {
    core.startGroup("Update article cache");
    exec(`python ${ACTION_PATH}/src/cache.py --update`);
    const article_cacheId = await cache.saveCache(article_paths, article_key);
    core.endGroup();
  } catch (error) {
    core.warning(error.message);
  }

  return output;
}
