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
var article_key = `article-${RUNNER_OS}-${ARTICLE_CACHE_NUMBER}-${randomId}`;
const article_restoreKeys = [`article-${RUNNER_OS}-${ARTICLE_CACHE_NUMBER}`];
const article_paths = [
  ".snakemake",
  ".showyourwork/tmp",
  ".showyourwork/resources",
  "environment.yml",
  "ms.pdf",
  "src",
];

// If this is one of the test repos, don't restore the cache!
const GITHUB_SLUG = shell.env["GITHUB_REPOSITORY"];
if (
  GITHUB_SLUG == "rodluger/showyourwork-template-minimal-test" ||
  GITHUB_SLUG == "rodluger/showyourwork-template-full-test"
) {
  article_key = `${randomId}`;
}

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
    exec("snakemake -c1 --use-conda --verbose --reason --notemp ms.pdf");
  } else {
    exec("snakemake -c1 --use-conda --reason --notemp ms.pdf");
  }
  output.push("ms.pdf");
  core.endGroup();

  // Build arxiv tarball
  if (core.getInput("arxiv-tarball") == "true") {
    const arxiv_tarball_exclude = getInputAsArray("arxiv-tarball-exclude").join(
      ","
    );
    core.startGroup("Build ArXiV tarball");
    if (core.getInput("verbose") == "true") {
      exec(
        `snakemake -c1 --use-conda --verbose --reason --notemp arxiv.tar.gz --config arxiv_tarball_exclude=${arxiv_tarball_exclude}`
      );
    } else {
      exec(
        `snakemake -c1 --use-conda --reason --notemp arxiv.tar.gz --config arxiv_tarball_exclude=${arxiv_tarball_exclude}`
      );
    }
    output.push("arxiv.tar.gz");
    core.endGroup();
  }

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
