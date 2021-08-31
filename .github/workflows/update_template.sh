#!/bin/bash

# Exit on error
set -e

# Update the template
git clone --recurse-submodules https://github.com/${USER}/${REPO}
cd ${REPO}/showyourwork
git fetch --all --tags
git checkout tags/"${{steps.bump.outputs.result}}"
cd ..
git add showyourwork
git -c user.name='gh-actions' -c user.email='gh-actions' commit -m "update showyourwork to ${{steps.bump.outputs.result}}"
git push https://x-access-token:${ACCESS_TOKEN}@github.com/${USER}/${REPO} main

# Now update the template-test repo; note that
# we sub in a different YAML file that will run the tests
cp showyourwork-action/resources/tests.yml .github/workflows/showyourwork.yml
git add .github/workflows/showyourwork.yml
git -c user.name='gh-actions' -c user.email='gh-actions' commit -m "update showyourwork to ${{steps.bump.outputs.result}}"
git push --force https://x-access-token:${ACCESS_TOKEN}@github.com/${USER}/${REPO}-test main
