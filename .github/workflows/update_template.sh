#!/bin/bash

# Exit on error
set -e

# Update the template
git clone --recurse-submodules --single-branch --branch ${BRANCH} https://github.com/rodluger/showyourwork-template
cd showyourwork-template/showyourwork
git fetch --all --tags
git checkout tags/${VERSION}
cd ..
git add showyourwork
git -c user.name='gh-actions' -c user.email='gh-actions' commit -m "update showyourwork to ${VERSION}"
git push https://x-access-token:${ACCESS_TOKEN}@github.com/showyourwork-template ${BRANCH}

# Now update the example repo (an instance of the template)
git push --force https://x-access-token:${ACCESS_TOKEN}@github.com/showyourwork-example ${BRANCH}
