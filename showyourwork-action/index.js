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
const CONDA_CACHE_NUMBER_DEV = 14;
const ARTICLE_CACHE_NUMBER_DEV = 12;
const CONDA_CACHE_NUMBER = core.getInput("conda-cache-number");
const ARTICLE_CACHE_NUMBER = core.getInput("article-cache-number");
const RUNNER_OS = shell.env["RUNNER_OS"];
const conda_key = `conda-${RUNNER_OS}-${CONDA_CACHE_NUMBER_DEV}-${CONDA_CACHE_NUMBER}`;
const conda_restoreKeys = [];
const conda_paths = ["~/.conda", "~/.condarc", "~/conda_pkgs_dir"];
const randomId = makeId(8);
const article_key = `article-${RUNNER_OS}-${ARTICLE_CACHE_NUMBER_DEV}-${ARTICLE_CACHE_NUMBER}-${randomId}`;
const article_restoreKeys = [
  `article-${RUNNER_OS}-${ARTICLE_CACHE_NUMBER_DEV}-${ARTICLE_CACHE_NUMBER}`,
];
const article_paths = [".showyourwork/cache", ".snakemake"];

// Repo we're running the action on
const GITHUB_WORKSPACE = shell.env["GITHUB_WORKSPACE"];
const GITHUB_REPOSITORY = shell.env["GITHUB_REPOSITORY"];
const GITHUB_USER = GITHUB_REPOSITORY.split("/")[0];
const CURRENT_BRANCH = shell
  .exec("git rev-parse --abbrev-ref HEAD")
  .replace(/(\r\n|\n|\r)/gm, "");
const GITHUB_TOKEN = core.getInput("github-token");

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
      const CONDA_URL = core.getInput("conda-url");
      exec(`wget --no-verbose ${CONDA_URL} -O ./conda.sh`, "Download conda");
      exec(
        "bash ./conda.sh -b -p ~/.conda && rm -f ./conda.sh",
        "Install conda"
      );
      exec(
        ". ~/.conda/etc/profile.d/conda.sh && " + 
        "conda config --add pkgs_dirs ~/conda_pkgs_dir"
      );
    }
    exec(". ~/.conda/etc/profile.d/conda.sh && conda info", "Conda info");
    
    // Create environment & install dependencies
    if (!shell.test("-d", "./envs")) {
      exec(
        ". ~/.conda/etc/profile.d/conda.sh && conda create -y -p ./envs",
        "Create environment"
      );
      exec_envs(
        "conda install -y -c defaults -c conda-forge -c bioconda " + 
        "mamba pip snakemake tectonic",
        "Install dependencies"
      );
    }

    // Install showyourwork from source
    var VERSION = core.getInput("showyourwork-version");
    if (VERSION == "latest") VERSION = "main";
    exec_envs(
      "python -m pip install " +
      `git+https://github.com/rodluger/showyourwork.git@${VERSION}`,
      "Install showyourwork"
    );
    core.endGroup();

    // Save conda cache
    try {
      core.startGroup("Update conda cache");
      const conda_cacheId = await cache.saveCache(conda_paths, conda_key);
      core.endGroup();
    } catch (error) {
      core.warning(error.message);
    }

    // Restore article cache
    core.startGroup("Restore article cache");
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
    core.startGroup("Build article");
    if (core.getInput("verbose") == "true") {
      exec_envs("showyourwork --verbose");
    } else {
      exec_envs("showyourwork");
    }
    outputs.push("ms.pdf");
    core.endGroup();

    // Generate DAG?
    if (core.getInput("generate-dag") == "true") {
      core.startGroup("Generate article DAG");
      exec_envs("showyourwork --dag");
      outputs.push("dag.pdf");
      core.endGroup();
    }

    // Save article cache
    try {
      core.startGroup("Update article cache");
      exec_envs("showyourwork --update-cache");
      const article_cacheId = await cache.saveCache(article_paths, article_key);
      core.endGroup();
    } catch (error) {
      core.warning(error.message);
    }

    // Upload artifact
    if (core.getInput("upload-artifact") == "true") {
      core.startGroup("Upload article artifact");
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
      core.startGroup("Uploading output");
      const TARGET_BRANCH = `${CURRENT_BRANCH}-pdf`;
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
        "git -c user.name='gh-actions' -c user.email='gh-actions' " + 
        "commit -m 'force-push article output'"
      );
      shell.exec(
        "git push --force " + 
        `https://x-access-token:${GITHUB_TOKEN}@github.com/${GITHUB_REPOSITORY} ` + 
        `${TARGET_BRANCH}`
      );
      core.endGroup();
    }

    // Format some files if this is a fresh repo based on the template
    if (shell.grep("<!--", "README.md").length > 1) {
      core.startGroup(`Format LICENSE and README.md`);
      shell.cd(GITHUB_WORKSPACE);
      // undo any current changes
      shell.exec("git reset --hard HEAD"); 
      // update in case things changed since the action started
      shell.exec(`git pull origin ${CURRENT_BRANCH}`); 
      shell.sed("-i", "Sit tight.*", "", "README.md");
      shell.sed("-i", "<!--", "", "README.md");
      shell.sed("-i", "-->", "", "README.md");
      shell.sed(
        "-i",
        "rodluger/showyourwork-template",
        GITHUB_REPOSITORY,
        "README.md"
      );
      shell.sed("-i", "Author Name", GITHUB_USER, "LICENSE");
      shell.exec("git add README.md");
      shell.exec("git add LICENSE");
      shell.exec(
        "git -c user.name='gh-actions' -c user.email='gh-actions' " + 
        "commit -m '[skip ci] format README.md'"
      );
      shell.exec(
        "git push " + 
        `https://x-access-token:${GITHUB_TOKEN}@github.com/${GITHUB_REPOSITORY} ` + 
        `${CURRENT_BRANCH}`
      );
      core.endGroup();
    }
  } catch (error) {
    core.setFailed(error.message);
  }
})();
