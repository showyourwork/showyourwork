#!/bin/bash

# Exit on error
set -e

# Update the template
git clone --recurse-submodules https://github.com/${USER}/${REPO}
cd ${REPO}/showyourwork
git fetch --all --tags
git checkout tags/${VERSION}
cd ..
git add showyourwork
git -c user.name='gh-actions' -c user.email='gh-actions' commit -m "update showyourwork to ${VERSION}"
git push https://x-access-token:${ACCESS_TOKEN}@github.com/${USER}/${REPO} main

# Now update the template-test repo; note that
# we sub in a different YAML file that will run the tests
cp ${GITHUB_WORKSPACE}/showyourwork-action/resources/tests.yml .github/workflows/showyourwork.yml
git add .github/workflows/showyourwork.yml
git -c user.name='gh-actions' -c user.email='gh-actions' commit -m "update showyourwork to ${VERSION}"
git push --force https://x-access-token:${ACCESS_TOKEN}@github.com/${USER}/${REPO}-test main
