#!/bin/bash
set -e

# Set defaults
rm -rf ~/.cookiecutterrc
echo "default_context:" >> ~/.cookiecutterrc
echo "    github_repo: \"${TARGET_REPOSITORY}\"" >> ~/.cookiecutterrc
echo "    github_user: \"${TARGET_USER}\"" >> ~/.cookiecutterrc
echo "    author_name: \"${AUTHOR_NAME}\"" >> ~/.cookiecutterrc
echo "    _showyourwork_sha: \"${SHOWYOURWORK_SHA}\"" >> ~/.cookiecutterrc
[[ ! -z "${ADD_CALLBACK}" ]] && echo "    _source_repo: \"${GITHUB_REPOSITORY}\"" >> ~/.cookiecutterrc

# Create repo and force push to target
cookiecutter --no-input ./showyourwork/cookiecutter-showyourwork
cd ${TARGET_REPOSITORY}
git init
git checkout --orphan main
git add .
git -c user.name="${GITHUB_REPOSITORY}" -c user.email="${GITHUB_REPOSITORY}" commit -m "auto commit from ${GITHUB_REPOSITORY}"
git push --force https://x-access-token:${ACCESS_TOKEN}@github.com/${TARGET_USER}/${TARGET_REPOSITORY} main
