# showyourwork-action

The **showyourwork-action** makes it easy to build open-source, reproducible scientific articles from a LaTeX manuscript and Python figure scripts. This action is intended to be run on repositories generated using the [showyourwork](https://github.com/rodluger/showyourwork) package.

This action generates all required figures from Python scripts in the `figures/` directory of your repo and compiles the LaTeX manuscript `tex/ms.tex` into the file `ms.pdf`, which is uploaded as an artifact upon completion of the workflow. This file is also automatically force-pushed to a branch on your repository with the same name as the current branch plus the suffix `-pdf`.

This action automatically caches all builds, so it won't re-run figure scripts that have not changed.

## Inputs:

### `action-path`

_Optional_. Path to this action relative to the top level of the repo. Default: `.showyourwork/showyourwork-action`

### `article-cache-number`

_Optional_. Bump this number to reset the article cache. Default: `0`

### `arxiv-tarball`

_Optional_. Build a tarball for easy ArXiV submission? This will be pushed to the same branch as the article output. Default `true`

### `conda-cache-number`

_Optional_. Bump this number to reset the `conda` cache. Default: `0`

### `conda-url`

_Optional_. Exact url pointing to the `conda` install script. Default: `https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh`

### `force-push`

_Optional_. Force-push output to branch `<xxx>-pdf`, where `<xxx>` is the current branch name? Default: `true`

### `generate-report`

_Optional_. Generate an HTML report of the build process and push it to GitHub Pages? Default: `true`

### `gh-pages-branch`

_Optional_. The branch serving GitHub Pages. Default: `gh-pages`

### `github-token`

_Optional_. A token for access to GitHub (e.g. `secrets.GITHUB_TOKEN`). Do not set this value explicitly -- always use a secret! Default: `${{ github.token }}` (usually set automatically)

### `verbose`

_Optional_. Enable verbose output and debug messages? Default: `false`

## Example usage

Below is a complete example that will automatically compile the article PDF on a repository with a `showyourwork`-conforming layout. Note that it is recommended to check out the repo with `fetch-depth: 0` to enable output caching in between commits.

```yaml
on: push
jobs:
  showyourwork:
    runs-on: ubuntu-latest
    name: Build the article PDF
    steps:
      - name: Checkout
        uses: actions/checkout@v2
        with:
          fetch-depth: 0
      - name: Build the article PDF
        id: build
        uses: ./.showyourwork/showyourwork-action
```
