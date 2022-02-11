// Imports
const core = require("@actions/core");
const cache = require("@actions/cache");
const shell = require("shelljs");
const constants = require("./constants.js");
const { makeId, exec } = require("./utils");

// Exports
module.exports = { buildArticle };

/**
 * Build the article.
 *
 */
async function buildArticle(ARTICLE_CACHE_NUMBER = null) {
  // Article cache settings. We only cache the contents of 
  // `.snakemake`, `.showyourwork`, and `src/tex/figures`.
  // Note that the GITHUB_REF (branch) is part of the cache key
  // so we don't mix up the caches for different branches!
  if (ARTICLE_CACHE_NUMBER == null)
    ARTICLE_CACHE_NUMBER = core.getInput("article-cache-number");
  const RUNNER_OS = shell.env["RUNNER_OS"];
  const GITHUB_REF = shell.env["GITHUB_REF"];
  const randomId = makeId(8);
  const article_key = `article-${constants.article_cache_version}-${RUNNER_OS}-${GITHUB_REF}-${ARTICLE_CACHE_NUMBER}-${randomId}`;
  const article_restoreKeys = [
    `article-${constants.article_cache_version}-${RUNNER_OS}-${GITHUB_REF}-${ARTICLE_CACHE_NUMBER}`,
  ];
  const article_paths = [
    ".showyourwork",
    ".snakemake",
    "src/tex/figures"
  ];

  // Restore the article cache
  core.startGroup("Restore article cache");
  const article_cacheKey = await cache.restoreCache(
    article_paths,
    article_key,
    article_restoreKeys
  );

  exec("make _restore_cache");
  core.endGroup();

  // Build the article
  core.startGroup("Build article");
  exec("make");
  core.endGroup();

  // Build arxiv tarball
  if (core.getInput("arxiv-tarball") == "true") {
    // TODO
  }

  // Save article cache
  core.startGroup("Update article cache");
  exec("make _update_cache");
  const article_cacheId = await cache.saveCache(article_paths, article_key);
  core.endGroup();
}
