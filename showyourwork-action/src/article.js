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
async function buildArticle(SHOWYOURWORK_VERSION, ARTICLE_CACHE_NUMBER = null) {
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
  const GITHUB_BRANCH = shell
    .exec("echo ${GITHUB_REF##*/}", {silent: true})
    .replace(/(\r\n|\n|\r)/gm, "");
  const GITHUB_SLUG = shell.env["GITHUB_REPOSITORY"];
  const randomId = makeId(8);
  const article_key = `article-${SHOWYOURWORK_VERSION}-${RUNNER_OS}-${GITHUB_REF}-${ARTICLE_CACHE_NUMBER}-${randomId}`;
  const article_restoreKeys = [
    `article-${SHOWYOURWORK_VERSION}-${RUNNER_OS}-${GITHUB_REF}-${ARTICLE_CACHE_NUMBER}`,
  ];
  const article_paths = [
    ".snakemake",
    ".showyourwork",
    ".last-commit",
    "environment.yml",
    "ms.pdf",
    "src",
  ];

  // Is this a unit test run?
  const UNIT_TEST = (
    (GITHUB_SLUG == "rodluger/showyourwork-example") || 
    (GITHUB_SLUG == "rodluger/showyourwork-example-dev")
  );

  // DEBUG
  core.info("ARTICLE_CACHE_NUMBER:");
  core.info(ARTICLE_CACHE_NUMBER);
  core.info("ARTICLE_CACHE_NUMBER == null:");
  core.info(ARTICLE_CACHE_NUMBER == null);
  core.info("ARTICLE_CACHE_NUMBER == '':");
  core.info(ARTICLE_CACHE_NUMBER == "");

  // We'll cache the article unless it's a unit test run
  // or if the user set the cache number to `null` (or empty).
  // But for good measure, we test the caching feature on 
  // the ``simple-figure`` branch when running unit tests
  const CACHE_ARTICLE = (
    !(
      UNIT_TEST || 
      ARTICLE_CACHE_NUMBER == null || 
      ARTICLE_CACHE_NUMBER == ""
    ) || 
    (
      GITHUB_BRANCH == "simple-figure"
    )
  );

  // Restore the article cache
  if (CACHE_ARTICLE) {
    core.startGroup("Restore article cache");
    const article_cacheKey = await cache.restoreCache(
      article_paths,
      article_key,
      article_restoreKeys
    );
    exec(`rm -f .showyourwork/repo.json`); // Always re-generate this!
    exec(`python ${ACTION_PATH}/src/cache.py --restore`);
    core.endGroup();
  }

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

  // Save article cache (failure OK)
  if (CACHE_ARTICLE) {
    try {
      core.startGroup("Update article cache");
      exec(`python ${ACTION_PATH}/src/cache.py --update`);
      const article_cacheId = await cache.saveCache(article_paths, article_key);
      core.endGroup();
    } catch (error) {
      core.warning(error.message);
    }
  }
  
  return output;
}
