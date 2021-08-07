// Imports
const core = require("@actions/core");
const artifact = require("@actions/artifact");
const github = require("@actions/github");
const cache = require("@actions/cache");
const shell = require("shelljs");

// Random hash
// https://stackoverflow.com/a/1349426
// https://github.com/actions/cache/issues/432#issuecomment-740376179
function makeId(length) {
  var result = "";
  var characters = "abcdefghijklmnopqrstuvwxyz0123456789";
  var charactersLength = characters.length;
  for (var i = 0; i < length; i++) {
    result += characters.charAt(Math.floor(Math.random() * charactersLength));
  }
  return result;
}

// Cache settings
const CONDA_CACHE_NUMBER_DEV = 13;
const ARTICLE_CACHE_NUMBER_DEV = 11;
const CONDA_CACHE_NUMBER = core.getInput("conda-cache-number");
const ARTICLE_CACHE_NUMBER = core.getInput("article-cache-number");
const RUNNER_OS = shell.env["RUNNER_OS"];
const conda_key = `conda-${RUNNER_OS}-${CONDA_CACHE_NUMBER_DEV}-${CONDA_CACHE_NUMBER}`;
const conda_restoreKeys = [];
const conda_paths = ["~/conda_pkgs_dir"];
const randomId = makeId(8);
const article_key = `article-${RUNNER_OS}-${ARTICLE_CACHE_NUMBER_DEV}-${ARTICLE_CACHE_NUMBER}-${randomId}`;
const article_restoreKeys = [
  `article-${RUNNER_OS}-${ARTICLE_CACHE_NUMBER_DEV}-${ARTICLE_CACHE_NUMBER}`,
];
const article_paths = [".showyourwork/cache", ".snakemake"];

// Exec, exit on failure
function exec(cmd, group) {
  if (typeof group !== "undefined") {
    core.startGroup(group);
  }
  if (shell.exec(cmd, { shell: "/bin/bash" }).code != 0) {
    shell.echo(`Error: ${cmd}`);
    shell.exit(1);
  }
  if (typeof group !== "undefined") {
    core.endGroup();
  }
}

// Exec in conda env `./envs`
function exec_envs(cmd, group) {
  exec(
    `. ~/.conda/etc/profile.d/conda.sh && conda activate ./envs && ${cmd}`,
    group
  );
}

(async () => {
  try {
    shell.set("-e");

    shell.echo("DEBUG!");
    shell.exec("ls");

    //
    if (shell.exec("grep '<!--' README.md", { shell: "/bin/bash" }).code == 0) {
      shell.echo("DEBUG!");
      core.startGroup(`Formatting README.md`);
      exec(`sed -i '' 's/<!--//' README.md`);
      exec(
        `sed -i '' 's/rodluger/showyourwork-template/${GITHUB_REPOSITORY}/' README.md`
      );
      exec(`sed -i '' 's/-->//' README.md`);
      exec(`git add README.md`);
      exec(
        "git -c user.name='gh-actions' -c user.email='gh-actions' commit -m 'format README.md'"
      );
      exec(`git push`);
      core.endGroup();
    }

    // Restore conda cache
    core.startGroup(`Restore conda cache`);
    const conda_cacheKey = await cache.restoreCache(
      conda_paths,
      conda_key,
      conda_restoreKeys
    );
    core.endGroup();

    // Setup conda
    if (!shell.test("-d", "~/.conda")) {
      const CONDA_URL = core.getInput("conda-url");
      exec(`wget --no-verbose ${CONDA_URL} -O ./conda.sh`, "Download conda");
      exec(
        "bash ./conda.sh -b -p ~/.conda && rm -f ./conda.sh",
        "Install conda"
      );
    }
    core.startGroup(`Conda info`);
    exec(
      ". ~/.conda/etc/profile.d/conda.sh && conda config --add pkgs_dirs ~/conda_pkgs_dir"
    );
    exec(". ~/.conda/etc/profile.d/conda.sh && conda info");
    core.endGroup();
    if (!shell.test("-d", "./envs")) {
      exec(
        ". ~/.conda/etc/profile.d/conda.sh && conda create -y -p ./envs",
        "Create environment"
      );
    }
    // Install showyourwork
    // TODO: If `showyourwork-version` is not a branch name or a SHA, install directly from conda-forge!
    var VERSION = core.getInput("showyourwork-version");
    if (VERSION == "latest") VERSION = "main";
    exec_envs(
      "conda install -y -c defaults -c conda-forge -c bioconda mamba pip snakemake tectonic",
      "Install dependencies"
    );
    exec_envs(
      `python -m pip install git+https://github.com/rodluger/showyourwork.git@${VERSION}`,
      "Install showyourwork"
    );
    core.endGroup();

    // Save conda cache
    try {
      core.startGroup(`Update conda cache`);
      const conda_cacheId = await cache.saveCache(conda_paths, conda_key);
      core.endGroup();
    } catch (error) {
      core.warning(error.message);
    }

    // Restore article cache
    core.startGroup(`Restore article cache`);
    const article_cacheKey = await cache.restoreCache(
      article_paths,
      article_key,
      article_restoreKeys
    );
    exec_envs("showyourwork --restore-cache");
    core.endGroup();

    // Outputs
    var outputs = [];

    // Build the article
    core.startGroup(`Build article`);
    if (core.getInput("verbose") == "true") {
      exec_envs("showyourwork --verbose");
    } else {
      exec_envs("showyourwork");
    }
    outputs.push("ms.pdf");
    core.endGroup();

    // Generate DAG?
    if (core.getInput("generate-dag") == "true") {
      core.startGroup(`Generate article DAG`);
      exec_envs("showyourwork --dag");
      outputs.push("dag.pdf");
      core.endGroup();
    }

    // Save article cache
    try {
      core.startGroup(`Update article cache`);
      exec_envs("showyourwork --update-cache");
      const article_cacheId = await cache.saveCache(article_paths, article_key);
      core.endGroup();
    } catch (error) {
      core.warning(error.message);
    }

    // Upload artifact
    if (core.getInput("upload-artifact") == "true") {
      core.startGroup(`Upload article artifact`);
      const artifactClient = artifact.create();
      const artifactName = "article-pdf";
      const rootDirectory = ".";
      const options = {
        continueOnError: false,
      };
      const uploadResponse = await artifactClient.uploadArtifact(
        artifactName,
        outputs,
        rootDirectory,
        options
      );
      core.endGroup();
    }

    // Force-push to `-pdf` branch
    if (core.getInput("force-push") == "true") {
      core.startGroup(`Uploading output`);
      const CURRENT_BRANCH = shell
        .exec("git rev-parse --abbrev-ref HEAD")
        .replace(/(\r\n|\n|\r)/gm, "");
      const TARGET_BRANCH = `${CURRENT_BRANCH}-pdf`;
      const GITHUB_TOKEN = core.getInput("github-token");
      const GITHUB_REPOSITORY = shell.env["GITHUB_REPOSITORY"];
      const TARGET_DIRECTORY = shell
        .exec("mktemp -d")
        .replace(/(\r\n|\n|\r)/gm, "");
      shell.cp("-R", ".", `${TARGET_DIRECTORY}`);
      shell.cd(`${TARGET_DIRECTORY}`);
      shell.exec(`git checkout --orphan ${TARGET_BRANCH}`);
      shell.exec("git rm --cached -rf .");
      for (const output of outputs) {
        shell.exec(`git add -f ${output}`);
      }
      shell.exec(
        "git -c user.name='gh-actions' -c user.email='gh-actions' commit -m 'force-push article output'"
      );
      shell.exec(
        `git push --force https://x-access-token:${GITHUB_TOKEN}@github.com/${GITHUB_REPOSITORY} ${TARGET_BRANCH}`
      );
      core.endGroup();
    }
  } catch (error) {
    core.setFailed(error.message);
  }
})();
